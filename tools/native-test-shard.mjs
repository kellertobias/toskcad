import {execFileSync} from 'node:child_process';
import {dirname,resolve} from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';

export const exampleGeneratorFilter='generator_rebuilds_the_shipped_project';
// Pinned cargo-nextest release; CI installs exactly this version.
export const nextestVersion='0.9.152';
// CI compiles the workspace tests once per platform into a nextest archive,
// then runs it on this many runners. `hash:K/N` partitions are disjoint and
// together select every test exactly once.
export const testPartitions=4;
export function archiveName(target){
  if(!/^[a-z0-9]+-[a-z0-9]+$/.test(target))throw new Error(`Invalid archive target: ${target}`);
  return `nextest-${target}.tar.zst`;
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
export function workspaceArchiveArguments(file){
  return ['nextest','archive','--workspace','--cargo-profile','ci','--locked','--archive-file',file];
}
// The archive is extracted so its `target/` lands where it was built (the
// workspace's own target directory). Tests that launch built executables
// use compile-time paths into that directory.
export function archiveShardArguments(file,partition,total=testPartitions){
  if(!Number.isInteger(total)||total<1)throw new Error(`Invalid partition count: ${total}`);
  if(!Number.isInteger(partition)||partition<1||partition>total)throw new Error(`Invalid test partition: ${partition}/${total}`);
  return ['nextest','run','--archive-file',file,'--workspace-remap','.','--extract-to','.','--extract-overwrite','--partition',`hash:${partition}/${total}`,'--no-fail-fast'];
}
// nextest does not run doctests; keep them in the suite with libtest.
export function workspaceDocTestArguments(){
  return ['test','--workspace','--doc','--profile','ci','--locked'];
}
// Local packaging tests the full workspace first. In CI the package job runs
// with `gated`: the archive shards and doctests of the same commit gate
// publication, so the package job builds and smoke-checks only.
export function packageTestArguments(mode='full'){
  if(mode==='full')return workspaceTestArguments();
  if(mode==='gated')return null;
  throw new Error(`Unknown package test mode: ${mode}`);
}
export function commandArguments(argv){
  const [mode='full',file,partition]=argv;
  if(mode==='archive'){
    if(!file)throw new Error('Usage: native-test-shard.mjs archive FILE');
    return workspaceArchiveArguments(file);
  }
  if(mode==='shard'){
    if(!file||!partition)throw new Error('Usage: native-test-shard.mjs shard FILE PARTITION');
    return archiveShardArguments(file,Number(partition));
  }
  if(file)throw new Error('Unexpected arguments');
  return workspaceTestArguments(mode);
}

if(process.argv[1]&&import.meta.url===pathToFileURL(resolve(process.argv[1])).href){
  try{
    execFileSync('cargo',commandArguments(process.argv.slice(2)),{
      cwd:resolve(dirname(fileURLToPath(import.meta.url)),'..'),stdio:'inherit',
    });
  }catch(error){
    console.error(error.message);
    process.exitCode=error.status??1;
  }
}
