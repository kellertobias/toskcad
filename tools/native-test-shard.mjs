import {execFileSync} from 'node:child_process';
import {mkdirSync,readFileSync,writeFileSync} from 'node:fs';
import {dirname,join,resolve} from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';

const root=resolve(dirname(fileURLToPath(import.meta.url)),'..');
export const exampleGeneratorFilter='generator_rebuilds_the_shipped_project';
// Pinned cargo-nextest release; CI installs exactly this version.
export const nextestVersion='0.9.152';
// CI compiles every workspace test once per platform, then splits the test
// binaries into this many archives, one per test runner. Every binary is in
// exactly one archive, so the runners together run every test exactly once.
export const testPartitions=4;
export function partitionArchiveName(partition,total=testPartitions){
  requirePartition(partition,total);
  return `nextest-${partition}.tar.zst`;
}
function requirePartition(partition,total){
  if(!Number.isInteger(total)||total<1)throw new Error(`Invalid partition count: ${total}`);
  if(!Number.isInteger(partition)||partition<1||partition>total)throw new Error(`Invalid test partition: ${partition}/${total}`);
}
// nextest runs every test in its own process across all cores. `test(...)`
// is a substring match, so the generator filter and its negation partition
// the suite exactly.
export function workspaceTestArguments(shard='full'){
  const args=['nextest','run','--workspace','--cargo-profile','ci','--locked','--no-fail-fast'];
  if(shard==='generators')args.push('-E',`test(${exampleGeneratorFilter})`);
  else if(shard==='remaining')args.push('-E',`not test(${exampleGeneratorFilter})`);
  else if(shard!=='full')throw new Error(`Unknown native test shard: ${shard}`);
  return args;
}
// Builds every workspace test and lists each binary with its tests.
export function workspaceListArguments(){
  return ['nextest','list','--workspace','--cargo-profile','ci','--locked','--message-format','json'];
}
export function partitionArchiveArguments(binaries,file){
  if(!binaries.length)throw new Error('An archive partition needs at least one test binary');
  const expression=binaries.map(id=>`binary_id(=${id})`).join(' | ');
  return ['nextest','archive','--workspace','--cargo-profile','ci','--locked','-E',expression,'--archive-file',file];
}
// The archive is extracted so its `target/` lands where it was built (the
// workspace's own target directory). Tests that launch built executables
// use compile-time paths into that directory.
export function archiveRunArguments(file){
  return ['nextest','run','--archive-file',file,'--workspace-remap','.','--extract-to','.','--extract-overwrite','--no-fail-fast'];
}
// Cargo creates CARGO_TARGET_TMPDIR (`target/tmp`) while building, but nextest
// archives do not contain it. Tests write scratch files and ledgers there, so
// a runner creates it before running an extracted archive.
export function prepareArchiveRun(workspace=root){
  const tmp=join(archiveTargetDirectory(workspace),'tmp');
  mkdirSync(tmp,{recursive:true});
  return tmp;
}
// Deterministic longest-first assignment: each binary goes to the partition
// with the least weight so far. Weights only balance runtime; every listed
// binary is assigned to exactly one partition.
export function partitionBinaries(suites,weights={},total=testPartitions){
  requirePartition(1,total);
  const ids=Object.keys(suites).sort();
  const weight=id=>weights[id]??Math.max(1,Object.keys(suites[id]?.testcases??{}).length);
  const order=[...ids].sort((a,b)=>weight(b)-weight(a)||(a<b?-1:a>b?1:0));
  const bins=Array.from({length:total},()=>({binaries:[],weight:0}));
  for(const id of order){
    let best=0;
    for(let k=1;k<total;k++)if(bins[k].weight<bins[best].weight)best=k;
    bins[best].binaries.push(id);bins[best].weight+=weight(id);
  }
  for(const bin of bins)bin.binaries.sort();
  const assigned=bins.flatMap(bin=>bin.binaries);
  if(assigned.length!==ids.length||new Set(assigned).size!==ids.length)throw new Error('Test binaries were not partitioned exactly once');
  return bins;
}
function cargo(args,options={}){
  return execFileSync('cargo',args,{cwd:root,maxBuffer:1<<30,...options});
}
// nextest stores the target directory as `target/`, and `archiveRunArguments`
// extracts it to the workspace's `target/`. Tests launch `multicad` through
// the compile-time CARGO_BIN_EXE path, so the archived binaries must be
// compiled in exactly that directory, not the configured `target/cargo`.
export function archiveTargetDirectory(workspace=root){
  return join(workspace,'target');
}
export function archivePartitions(directory,total=testPartitions){
  const targetDirectory=process.env.CARGO_TARGET_DIR&&resolve(process.env.CARGO_TARGET_DIR);
  const same=(a,b)=>process.platform==='win32'?a?.toLowerCase()===b.toLowerCase():a===b;
  if(!same(targetDirectory,archiveTargetDirectory()))throw new Error(`Set CARGO_TARGET_DIR=${archiveTargetDirectory()} to archive tests: extracted archives place binaries there`);
  const listed=JSON.parse(cargo(workspaceListArguments(),{encoding:'utf8',stdio:['ignore','pipe','inherit']}));
  const suites=listed['rust-suites'];
  const weights=JSON.parse(readFileSync(join(root,'tools/test-weights.json'),'utf8')).seconds;
  const bins=partitionBinaries(suites,weights,total);
  mkdirSync(directory,{recursive:true});
  const manifest={partitions:bins.map((bin,i)=>({archive:partitionArchiveName(i+1,total),weight:Math.round(bin.weight),binaries:bin.binaries}))};
  for(const part of manifest.partitions){
    console.log(`${part.archive}: ${part.binaries.length} binaries, weight ${part.weight}`);
    cargo(partitionArchiveArguments(part.binaries,join(directory,part.archive)),{stdio:'inherit'});
  }
  writeFileSync(join(directory,'partitions.json'),JSON.stringify(manifest,null,1)+'\n');
  return manifest;
}
// nextest does not run doctests; keep them in the suite with libtest.
export function workspaceDocTestArguments(){
  return ['test','--workspace','--doc','--profile','ci','--locked'];
}
// Local packaging tests the full workspace first. In CI the package job runs
// with `gated`: the archive runs and doctests of the same commit gate
// publication, so the package job builds and smoke-checks only.
export function packageTestArguments(mode='full'){
  if(mode==='full')return workspaceTestArguments();
  if(mode==='gated')return null;
  throw new Error(`Unknown package test mode: ${mode}`);
}

if(process.argv[1]&&import.meta.url===pathToFileURL(resolve(process.argv[1])).href){
  const [mode='full',path,...extra]=process.argv.slice(2);
  try{
    if(extra.length)throw new Error('Unexpected arguments');
    if(mode==='archive'){
      if(!path)throw new Error('Usage: native-test-shard.mjs archive DIRECTORY');
      archivePartitions(resolve(path));
    }else if(mode==='shard'){
      if(!path)throw new Error('Usage: native-test-shard.mjs shard ARCHIVE');
      prepareArchiveRun();
      cargo(archiveRunArguments(path),{stdio:'inherit'});
    }else{
      if(path)throw new Error('Unexpected arguments');
      cargo(workspaceTestArguments(mode),{stdio:'inherit'});
    }
  }catch(error){
    console.error(error.message);
    process.exitCode=error.status??1;
  }
}
