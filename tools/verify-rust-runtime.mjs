// Developer and CI check that a compiled Rust-only executable meets the
// production runtime contract: Rust geometry identity, no OpenCascade or
// TK library linkage, no statically linked OpenCascade code, and only
// system runtime libraries. Usage:
//   node tools/verify-rust-runtime.mjs target/cargo/release/multicad
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
import {resolve} from 'node:path';
import {fileURLToPath} from 'node:url';
import {auditRuntime} from './native-package.mjs';

// Mirrors the archive verifier in package-release.mjs.
const occtLibrary=/(?:^|[/\\\s])(?:lib)?(?:tk[a-z0-9]|occt|opencascade)/i;
// Class names that any linked OpenCascade modeling, STEP or meshing code
// carries in its C++ runtime type information, even in stripped binaries.
export const occtMarkers=Object.freeze(['Standard_Transient','TopoDS_TShape','BRepCheck_Analyzer','STEPControl_Reader','BRepMesh_IncrementalMesh','Geom_BSplineSurface']);

export function checkDiagnostics(text){
  const lines=text.split(/\r?\n/);
  const kernel=lines.filter(line=>line.startsWith('CAD kernel: '));
  if(kernel.length!==1||kernel[0]!=='CAD kernel: MultiCAD Rust')throw new Error(`Executable does not identify the Rust kernel: ${kernel.join(' | ')||'missing'}`);
  if(lines.some(line=>line.startsWith('OCCT source:')||/OCCT \d/.test(line)))throw new Error('Executable diagnostics report OpenCascade');
  if(!lines.some(line=>line.startsWith('Runtime: compiled Rust')))throw new Error('Executable diagnostics do not state the compiled Rust runtime');
}
export function checkLibraries(libraries){
  if(!Array.isArray(libraries)||libraries.length===0)throw new Error('No runtime libraries were reported');
  const linked=libraries.filter(library=>occtLibrary.test(library));
  if(linked.length)throw new Error(`OpenCascade runtime libraries are linked: ${linked.join(', ')}`);
}
export function checkStaticCode(bytes){
  const found=occtMarkers.filter(marker=>bytes.includes(Buffer.from(marker)));
  if(found.length)throw new Error(`OpenCascade code is statically linked: ${found.join(', ')}`);
}
export function verifyRustRuntime(binary){
  const executable=resolve(binary);
  const diagnostics=execFileSync(executable,['--diagnostics'],{encoding:'utf8'});
  checkDiagnostics(diagnostics);
  const libraries=auditRuntime(executable);
  checkLibraries(libraries);
  checkStaticCode(readFileSync(executable));
  return {executable,libraries,kernel:'MultiCAD Rust'};
}
if(process.argv[1]&&resolve(process.argv[1])===fileURLToPath(import.meta.url)){
  if(process.argv.length!==3){console.error('Usage: node tools/verify-rust-runtime.mjs EXECUTABLE');process.exit(2);}
  console.log(JSON.stringify(verifyRustRuntime(process.argv[2]),null,2));
}
