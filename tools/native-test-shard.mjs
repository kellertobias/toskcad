import {execFileSync} from 'node:child_process';
import {dirname,resolve} from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';

export const exampleGeneratorFilter='generator_rebuilds_the_shipped_project';
// Pinned cargo-nextest release; CI installs exactly this version.
export const nextestVersion='0.9.152';
// nextest runs every test in its own process across all cores. `test(...)`
// is a substring match, so the generator and remaining shards partition the
// suite exactly as the former libtest filter and `--skip` pair did.
export function workspaceTestArguments(shard='full'){
  const args=['nextest','run','--workspace','--cargo-profile','ci','--locked','--no-fail-fast'];
  if(shard==='generators')args.push('-E',`test(${exampleGeneratorFilter})`);
  else if(shard==='remaining')args.push('-E',`not test(${exampleGeneratorFilter})`);
  else if(shard!=='full')throw new Error(`Unknown native test shard: ${shard}`);
  return args;
}
// nextest does not run doctests; keep them in the suite with libtest.
export function workspaceDocTestArguments(){
  return ['test','--workspace','--doc','--profile','ci','--locked'];
}
export function packageTestArguments(shard='full',platform=process.platform){
  if(shard!=='full'&&!(platform==='win32'&&shard==='remaining'))throw new Error('Only the paired Windows package job may select the remaining test shard');
  return workspaceTestArguments(shard);
}

if(process.argv[1]&&import.meta.url===pathToFileURL(resolve(process.argv[1])).href){
  try{
    execFileSync('cargo',workspaceTestArguments(process.argv[2]??'full'),{
      cwd:resolve(dirname(fileURLToPath(import.meta.url)),'..'),stdio:'inherit',
    });
  }catch(error){
    console.error(error.message);
    process.exitCode=error.status??1;
  }
}
