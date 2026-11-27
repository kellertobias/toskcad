import test from 'node:test';
import assert from 'node:assert/strict';
import {summarizeResize} from './measure-window-resize.mjs';

const trace = begins => [
  {phase:'start',version:1,t:0.5,unix_time:1000.5},
  ...begins.map((t,index)=>({phase:'frame-begin',frame:index+1,t})),
].map(row=>JSON.stringify(row)).join('\n');
// Epochs 1001 and 1002 are trace times 1.0 s and 2.0 s.
const driver = 'window x=0 y=0 width=1360 height=860\nresize-start 1001.000000\nresize-end 1002.000000\n';

test('keeps only frames that began during the resize', () => {
  const summary = summarizeResize(trace([0.2, 1.0, 1.01, 1.02, 1.04, 2.5]), driver);
  assert.equal(summary.frames, 4);
  assert.equal(summary.interval.count, 3);
  assert.ok(Math.abs(summary.interval.max - 20) < 1e-6);
  assert.equal(summary.ok, true);
});

test('fails on a drawable-timeout stall during the resize', () => {
  const summary = summarizeResize(trace([1.0, 1.008, 1.016, 1.9, 1.91]), driver);
  assert.equal(summary.slow, 1);
  assert.equal(summary.ok, false);
  assert.equal(summarizeResize(trace([1.0, 1.008, 1.016, 1.9, 1.91]), driver, 1000).ok, true);
});

test('rejects traces and driver output it cannot place', () => {
  assert.throws(() => summarizeResize('{"phase":"frame-begin","t":1}', driver), /start record/);
  assert.throws(() => summarizeResize(trace([1.1, 1.2]), 'resize-start 1001\n'), /resize-end/);
  assert.throws(() => summarizeResize(trace([1.1]), driver), /Fewer than two frames/);
  assert.throws(() => summarizeResize(trace([1.1, 1.2]), 'resize-start 1002\nresize-end 1001\n'), /ends before/);
});
