// Developer-only reader for MULTICAD_RESPONSE_TRACE files (CAD-347). The app
// writes them only when that variable is set; see
// crates/multicad-core/src/response_trace.rs for the record contract.
//
//   node tools/measure-response-trace.mjs TRACE.jsonl [--from SECONDS] [--through SECONDS]
//
// Times are in-process: an input starts when a frame first sees it and a
// result is visible when that frame's commands are encoded for presentation.
// Operating-system queueing before the frame and display scan-out after
// encoding are outside this trace; measure-native-input.mjs covers them.
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';

const STAGES = ['preparation','compilation','query','publication','upload'];
const ACTIONS = new Set(['pointer-down','key-down','text','paste','scroll','zoom','accessibility-action']);
const ms = seconds => seconds*1000;
const distribution = values => {
  if (!values.length) return null;
  const sorted = [...values].sort((a,b)=>a-b);
  const at = q => sorted[Math.ceil(sorted.length*q)-1];
  return {count:sorted.length,p50:at(.5),p95:at(.95),max:sorted.at(-1)};
};

export function parseResponseTrace(text) {
  const rows = text.trim().split('\n').filter(Boolean).map((line,index)=>{
    try { return JSON.parse(line); } catch { throw new Error(`Invalid JSON at line ${index+1}`); }
  });
  if (rows[0]?.phase!=='start' || rows[0].version!==1) throw new Error('Missing or unsupported response trace start record');
  for (const row of rows) {
    const times = row.phase==='span' ? [row.start,row.end] : [row.t];
    if (!times.every(Number.isFinite) || (row.phase==='span' && row.end<row.start)) throw new Error(`Invalid time in ${row.phase} record`);
    if (row.phase==='span' && !STAGES.includes(row.stage)) throw new Error(`Unknown stage ${row.stage}`);
  }
  return rows;
}

function frameTable(rows) {
  const frames = new Map();
  const get = frame => frames.get(frame) ?? frames.set(frame,{frame}).get(frame);
  for (const row of rows) {
    if (row.phase==='frame-begin') get(row.frame).begin=row.t;
    if (row.phase==='frame-built') get(row.frame).built=row.t;
    if (row.phase==='frame-encoded') get(row.frame).encoded ??= row.t;
  }
  return [...frames.values()].filter(f=>Number.isFinite(f.begin)&&Number.isFinite(f.built)).sort((a,b)=>a.begin-b.begin);
}

// The frame whose build contains `t`, else the first frame that begins after it.
function frameAt(frames, t) {
  return frames.find(f=>f.begin<=t && t<=f.built) ?? frames.find(f=>f.begin>t);
}
// When a result handled at `t` reaches an encoded frame.
function visibleAfter(frames, t) {
  const index = frames.findIndex(f=>f.begin<=t && t<=f.built);
  const candidates = index<0 ? frames.filter(f=>f.begin>t) : frames.slice(index);
  return candidates.find(f=>Number.isFinite(f.encoded));
}

const completion = {
  'scene-request': (rows,request)=>rows.find(r=>r.phase==='span' && r.label==='scene-accept' && r.generation===request.generation && r.start>=request.t),
  'model-preview-request': (rows,request)=>rows.find(r=>r.phase==='span' && r.label==='model-preview-accept' && r.generation===request.generation && r.start>=request.t),
  'drawing-request': (rows,request)=>rows.find(r=>r.phase==='drawing-reply' && r.revision===request.revision && r.t>=request.t),
  'preparation-request': (rows,request)=>rows.find(r=>r.phase==='preparation-reply' && r.generation===request.generation && r.t>=request.t),
};
const doneAt = row => row.phase==='span' ? row.end : row.t;
const succeeded = row => row.phase==='span' ? row.ok!==false : row.ok!==false && row.current!==false;

export function summarizeResponseTrace(text, {from=-Infinity, through=Infinity}={}) {
  const rows = parseResponseTrace(text);
  const frames = frameTable(rows);
  const inputs = rows.filter(r=>r.phase==='input' && ACTIONS.has(r.kind) && r.t>=from && r.t<=through);
  const feedback = [];
  const unpainted = [];
  for (const input of inputs) {
    const frame = frames.find(f=>f.frame===input.frame);
    if (!frame || !Number.isFinite(frame.encoded)) { unpainted.push(input.sequence); continue; }
    feedback.push(ms(frame.encoded-frame.begin));
  }
  const selectedFrames = frames.filter(f=>f.begin>=from && f.begin<=through);
  // A preparation reply that immediately requests a native preview is one
  // operator action: report the chain from the original request.
  const chains = [];
  const consumed = new Set();
  for (const request of rows.filter(r=>completion[r.phase] && r.t>=from && r.t<=through)) {
    if (consumed.has(request)) continue;
    let done = completion[request.phase](rows,request);
    const steps = [request.phase];
    while (done?.phase==='preparation-reply' && succeeded(done)) {
      const frame = frameAt(frames, done.t);
      const next = frame && rows.find(r=>r.phase==='model-preview-request' && r.t>=done.t && r.t<=frame.built);
      if (!next) break;
      consumed.add(next);
      steps.push(next.phase);
      done = completion[next.phase](rows,next);
    }
    const trigger = [...rows].reverse().find(r=>r.phase==='input' && ACTIONS.has(r.kind) && r.t<=request.t);
    const start = trigger ? frames.find(f=>f.frame===trigger.frame)?.begin ?? trigger.t : request.t;
    const chain = {request:steps.join(' > '), revision:request.revision ?? null, generation:request.generation ?? null,
      trigger:trigger ? {sequence:trigger.sequence, kind:trigger.kind, detail:trigger.detail ?? null, msBeforeRequest:ms(request.t-trigger.t)} : null};
    if (!done) { chains.push({...chain, outcome:'unanswered'}); continue; }
    const visible = visibleAfter(frames, doneAt(done));
    const end = visible?.encoded ?? doneAt(done);
    const stages = Object.fromEntries(STAGES.map(stage=>[stage,0]));
    const labels = {};
    for (const span of rows.filter(r=>r.phase==='span' && r.end>start && r.start<end)) {
      const overlap = ms(Math.min(span.end,end)-Math.max(span.start,start));
      const key = `${span.stage}:${span.label}@${span.thread ?? 'unnamed'}`;
      labels[key] = (labels[key] ?? 0)+overlap;
    }
    // Stage totals count each thread's outermost spans once, so a command span
    // that encloses its own scene request is not added twice.
    const spans = rows.filter(r=>r.phase==='span' && r.end>start && r.start<end);
    for (const span of spans) {
      const enclosed = spans.some(other=>other!==span && other.thread===span.thread && other.start<=span.start && other.end>=span.end && (other.start<span.start || other.end>span.end));
      if (!enclosed) stages[span.stage] += ms(Math.min(span.end,end)-Math.max(span.start,start));
    }
    const framesBetween = frames.filter(f=>f.begin>=start && f.built<=end);
    chains.push({...chain, outcome:succeeded(done) ? (visible ? 'visible' : 'handled, not painted') : 'failed',
      totalMs:ms(end-start), resultHandledMs:ms(doneAt(done)-start),
      paintMs: visible ? ms(visible.encoded-visible.built) : null,
      longestFrameBuildMs: framesBetween.length ? Math.max(...framesBetween.map(f=>ms(f.built-f.begin))) : null,
      stagesMs:stages, spansMs:labels});
  }
  const statuses = rows.filter(r=>r.phase==='frame-built' && r.status!==undefined && r.t>=from && r.t<=through).map(r=>({t:r.t, mode:r.mode, status:r.status}));
  const cancels = rows.filter(r=>['preparation-cancel','superseded'].includes(r.phase) && r.t>=from && r.t<=through).map(r=>({phase:r.phase, t:r.t}));
  return {
    scope:'in-process: frame that first sees an input to the encoded frame that shows its result; excludes OS queueing and display scan-out',
    build:{version:rows[0].product_version, commit:rows[0].commit ?? null, pid:rows[0].pid},
    frames:{count:selectedFrames.length,
      buildMs:distribution(selectedFrames.map(f=>ms(f.built-f.begin))),
      paintMs:distribution(selectedFrames.filter(f=>Number.isFinite(f.encoded)).map(f=>ms(f.encoded-f.built)))},
    inputs:{count:inputs.length, unpaintedSequences:unpainted, feedbackMs:distribution(feedback)},
    chains, statuses, cancels,
  };
}

if (import.meta.url===pathToFileURL(process.argv[1] ?? '').href) {
  const [path,...flags] = process.argv.slice(2);
  if (!path) { console.error('Usage: node tools/measure-response-trace.mjs TRACE.jsonl [--from SECONDS] [--through SECONDS]'); process.exit(2); }
  const option = name => { const index=flags.indexOf(name); return index<0 ? undefined : Number(flags[index+1]); };
  console.log(JSON.stringify(summarizeResponseTrace(readFileSync(path,'utf8'),{from:option('--from') ?? -Infinity, through:option('--through') ?? Infinity}),null,2));
}
