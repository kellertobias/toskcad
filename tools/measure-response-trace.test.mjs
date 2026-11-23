import assert from 'node:assert/strict';
import test from 'node:test';
import {summarizeResponseTrace} from './measure-response-trace.mjs';

const start = {phase:'start',t:0,version:1,pid:1,product_version:'0.1.0',commit:null};
const frame = (n,begin,built,encoded,extra={})=>[
  {phase:'frame-begin',t:begin,frame:n},
  {phase:'frame-built',t:built,frame:n,...extra},
  ...(encoded===undefined ? [] : [{phase:'frame-encoded',t:encoded,frame:n}]),
];
const span = (stage,label,s,e,fields={})=>({phase:'span',stage,label,start:s,end:e,thread:'main',...fields});
const trace = rows=>[start,...rows].map(row=>JSON.stringify(row)).join('\n');

test('an Apply input is measured through worker stages to the frame that shows the scene',()=>{
  const text = trace([
    ...frame(1,1.000,1.010,1.012),
    {phase:'input',t:1.000,sequence:1,frame:1,kind:'pointer-down'},
    span('preparation','command',1.002,1.008,{command:'Apply'}),
    {phase:'scene-request',t:1.007,revision:4,generation:9},
    span('compilation','model',1.010,1.200,{generation:9,thread:'multicad-model'}),
    span('publication','publish',1.200,1.230,{generation:9,thread:'multicad-model'}),
    ...frame(2,1.240,1.260,1.275,{status:'Parameters applied'}),
    span('publication','scene-accept',1.241,1.250,{generation:9,ok:true,published:true}),
    span('upload','viewport',1.262,1.270,{generation:9}),
  ]);
  const report = summarizeResponseTrace(text);
  assert.equal(report.chains.length,1);
  const [chain] = report.chains;
  assert.equal(chain.outcome,'visible');
  assert.equal(chain.trigger.sequence,1);
  assert.ok(Math.abs(chain.totalMs-275)<1e-6);
  assert.ok(Math.abs(chain.resultHandledMs-250)<1e-6);
  assert.ok(Math.abs(chain.stagesMs.compilation-190)<1e-6);
  // The command span on the main thread is counted once with its contents.
  assert.ok(Math.abs(chain.stagesMs.preparation-6)<1e-6);
  assert.ok(Math.abs(chain.stagesMs.publication-39)<1e-6);
  assert.ok(Math.abs(chain.stagesMs.upload-8)<1e-6);
  assert.ok(Math.abs(chain.paintMs-15)<1e-6);
  assert.deepEqual(report.statuses.map(s=>s.status),['Parameters applied']);
  assert.ok(Math.abs(report.inputs.feedbackMs.max-12)<1e-6);
});

test('a preparation that requests a native preview is one chain and a failure is not visible',()=>{
  const text = trace([
    ...frame(1,1,1.01,1.02),
    {phase:'input',t:1,sequence:1,frame:1,kind:'pointer-down'},
    {phase:'preparation-request',t:1.005,generation:3},
    span('preparation','model-request',1.006,1.3,{generation:3,thread:'model-preparation'}),
    ...frame(2,1.31,1.33,1.34),
    {phase:'preparation-reply',t:1.315,generation:3,current:true,ok:true},
    {phase:'model-preview-request',t:1.32,revision:0,generation:7},
    ...frame(3,1.9,1.95,1.96),
    span('publication','model-preview-accept',1.91,1.94,{generation:7,ok:false}),
  ]);
  const [chain] = summarizeResponseTrace(text).chains;
  assert.equal(chain.request,'preparation-request > model-preview-request');
  assert.equal(chain.outcome,'failed');
  assert.ok(Math.abs(chain.totalMs-960)<1e-6);
});

test('unanswered requests, unpainted inputs and cancellation are reported, not hidden',()=>{
  const text = trace([
    ...frame(1,1,1.01),
    {phase:'input',t:1,sequence:1,frame:1,kind:'key-down',detail:'Escape'},
    {phase:'drawing-request',t:1.005,revision:2,page:true},
    {phase:'preparation-cancel',t:1.006},
    {phase:'input',t:1.007,sequence:2,frame:1,kind:'key-up'},
  ]);
  const report = summarizeResponseTrace(text);
  assert.equal(report.chains[0].outcome,'unanswered');
  assert.deepEqual(report.inputs.unpaintedSequences,[1]);
  assert.equal(report.inputs.count,1);
  assert.deepEqual(report.cancels.map(c=>c.phase),['preparation-cancel']);
});

test('malformed traces are rejected',()=>{
  assert.throws(()=>summarizeResponseTrace('{"phase":"input","t":1}'),/start record/);
  assert.throws(()=>summarizeResponseTrace(trace([span('compilation','model',2,1)])),/Invalid time/);
  assert.throws(()=>summarizeResponseTrace(trace([span('paint','x',1,2)])),/Unknown stage/);
  assert.throws(()=>summarizeResponseTrace(`${JSON.stringify(start)}\n{oops`),/Invalid JSON at line 2/);
});

test('a time window selects inputs, frames and requests',()=>{
  const text = trace([
    ...frame(1,1,1.01,1.02),
    {phase:'input',t:1,sequence:1,frame:1,kind:'pointer-down'},
    {phase:'drawing-request',t:1.005,revision:1},
    ...frame(2,2,2.01,2.02),
    {phase:'drawing-reply',t:2.005,revision:1,ok:true,views:2,ready:2},
    ...frame(3,5,5.01,5.02),
    {phase:'input',t:5,sequence:2,frame:3,kind:'pointer-down'},
  ]);
  const early = summarizeResponseTrace(text,{through:3});
  assert.equal(early.inputs.count,1);
  assert.equal(early.frames.count,2);
  assert.equal(early.chains[0].outcome,'visible');
  assert.ok(Math.abs(early.chains[0].totalMs-1020)<1e-6);
  assert.equal(summarizeResponseTrace(text,{from:4}).chains.length,0);
});
