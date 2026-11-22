import assert from 'node:assert/strict';
import test from 'node:test';
import {checkDiagnostics,checkLibraries,checkStaticCode,occtMarkers} from './verify-rust-runtime.mjs';

const rust='ToskCAD 0.1.0\nStage: native-integration\nCAD kernel: MultiCAD Rust\nRuntime: compiled Rust and native libraries; no Python or Node.js dependency';
test('Rust diagnostics are accepted and OpenCascade identities are refused',()=>{
  checkDiagnostics(rust);
  assert.throws(()=>checkDiagnostics(rust.replace('MultiCAD Rust','OpenCascade 7.9.3')),/Rust kernel/);
  assert.throws(()=>checkDiagnostics(rust+'\nOCCT source: abc'),/OpenCascade/);
  assert.throws(()=>checkDiagnostics(rust.replace(/Runtime:.*/,'')),/compiled Rust runtime/);
});
test('dynamic OpenCascade and TK libraries are refused',()=>{
  checkLibraries(['/usr/lib/libSystem.B.dylib','/System/Library/Frameworks/Metal.framework/Versions/A/Metal']);
  assert.throws(()=>checkLibraries([]),/No runtime libraries/);
  for(const library of ['/usr/lib/libTKernel.dylib','/opt/occt/lib/libTKSTEP.7.9.dylib','C:\\occt\\TKBRep.dll','libopencascade.so'])
    assert.throws(()=>checkLibraries(['/usr/lib/libSystem.B.dylib',library]),/OpenCascade/);
});
test('statically linked OpenCascade type information is refused',()=>{
  checkStaticCode(Buffer.from('plain Rust executable bytes'));
  for(const marker of occtMarkers)
    assert.throws(()=>checkStaticCode(Buffer.from(`\0prefix${marker}suffix\0`)),/statically linked/);
});
