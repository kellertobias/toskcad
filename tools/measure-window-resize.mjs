// Developer-only window-resize frame probe (CAD-420). Node.js is not a
// runtime dependency. It reads a MULTICAD_RESPONSE_TRACE file and the output
// of tools/measure-window-resize.swift, keeps the frames that began during
// the resize, and fails when any interval between frame starts reaches the
// limit. Before CAD-420, frames painted during a resize could wait the full
// 1 s CAMetalLayer nextDrawable timeout.
//
//   HOME=$ISOLATED MULTICAD_RESPONSE_TRACE=$TMP/resize.jsonl target/release/multicad &
//   "$TMP/resize-driver" $! grip > $TMP/resize.log
//   node tools/measure-window-resize.mjs $TMP/resize.jsonl $TMP/resize.log [--max-interval MS]
//
// Times come from the application's frame-begin records, so the interval
// includes presentation waits. Place the window on the display under test
// first; refresh rate changes the expected median (8.3 ms at 120 Hz).
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';

const distribution = values => {
  if (!values.length) return null;
  const sorted = [...values].sort((a,b)=>a-b);
  const at = q => sorted[Math.min(sorted.length-1, Math.ceil(sorted.length*q)-1)];
  return {count:sorted.length,p50:at(.5),p95:at(.95),p99:at(.99),max:sorted.at(-1)};
};

export function summarizeResize(traceText, driverText, maxIntervalMs = 250) {
  const rows = traceText.trim().split('\n').filter(Boolean).map((line,index)=>{
    try { return JSON.parse(line); } catch { throw new Error(`Invalid JSON at trace line ${index+1}`); }
  });
  const start = rows[0];
  if (start?.phase!=='start' || start.version!==1 || !Number.isFinite(start.unix_time) || !Number.isFinite(start.t)) throw new Error('Missing or unsupported response trace start record');
  const epoch = label => {
    const match = driverText.match(new RegExp(`^${label} (\\d+(?:\\.\\d+)?)$`,'m'));
    if (!match) throw new Error(`Driver output has no ${label} time`);
    return Number(match[1]) - (start.unix_time - start.t);
  };
  const from = epoch('resize-start'), through = epoch('resize-end');
  if (!(through > from)) throw new Error('Resize ends before it starts');
  const begins = rows.filter(row=>row.phase==='frame-begin' && row.t>=from && row.t<=through).map(row=>row.t).sort((a,b)=>a-b);
  if (begins.length < 2) throw new Error('Fewer than two frames began during the resize');
  const intervals = begins.slice(1).map((t,index)=>(t-begins[index])*1000);
  const slow = intervals.filter(value=>value>=maxIntervalMs).length;
  return {seconds:through-from, frames:begins.length, interval:distribution(intervals), slow, maxIntervalMs, ok:slow===0};
}

if (import.meta.url===pathToFileURL(process.argv[1]??'').href) {
  const [trace, driver, ...options] = process.argv.slice(2);
  let limit = 250;
  if (options[0]==='--max-interval' && Number(options[1])>0) limit = Number(options[1]);
  else if (options.length) { console.error('Usage: measure-window-resize.mjs TRACE.jsonl DRIVER.log [--max-interval MS]'); process.exit(2); }
  if (!trace || !driver) { console.error('Usage: measure-window-resize.mjs TRACE.jsonl DRIVER.log [--max-interval MS]'); process.exit(2); }
  const summary = summarizeResize(readFileSync(trace,'utf8'), readFileSync(driver,'utf8'), limit);
  const ms = value => `${value.toFixed(1)} ms`;
  const {interval} = summary;
  console.log(`resize ${summary.seconds.toFixed(2)} s, ${summary.frames} frames`);
  console.log(`frame interval p50 ${ms(interval.p50)} p95 ${ms(interval.p95)} p99 ${ms(interval.p99)} max ${ms(interval.max)}`);
  console.log(`${summary.slow} interval(s) at or above ${limit} ms`);
  process.exit(summary.ok ? 0 : 1);
}
