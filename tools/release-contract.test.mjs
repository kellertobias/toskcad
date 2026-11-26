import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { createRequire } from 'node:module';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';
import { githubRepositoryUrl, requirePublicationEnvironment, requireVersion } from './release-policy.mjs';
import {archiveRunArguments,exampleGeneratorFilter,nextestVersion,packageTestArguments,partitionArchiveArguments,partitionArchiveName,partitionBinaries,testPartitions,workspaceDocTestArguments,workspaceListArguments,workspaceTestArguments} from './native-test-shard.mjs';

const require = createRequire(import.meta.url);
const config = require('../release.config.cjs');
const git = (cwd, args) => execFileSync('git', args, {cwd, encoding:'utf8', stdio:['ignore','pipe','pipe']}).trim();
function preview(sourceRoot) {
  // semantic-release hooks stdout; isolate it from node:test's worker protocol.
  const environment={...process.env};delete environment.NODE_TEST_CONTEXT;
  // The fixture is its own repository on main. Hide the runner's CI identity
  // so a test run on a feature branch does not read that branch as the target.
  for(const key of Object.keys(environment))if(key==='CI'||key.startsWith('GITHUB_'))delete environment[key];
  const program=`import { Writable } from 'node:stream';
    import { runRelease } from ${JSON.stringify(new URL('./run-semantic-release.mjs',import.meta.url).href)};
    const sink=new Writable({write(_chunk,_encoding,done){done();}});
    const result=await runRelease({dryRun:true,sourceRoot:process.argv[1],stdout:sink,stderr:sink});
    console.log(JSON.stringify(result===false?{version:null}:{version:result.nextRelease.version,notes:result.nextRelease.notes}));`;
  return JSON.parse(execFileSync(process.execPath,['--input-type=module','--eval',program,sourceRoot],{encoding:'utf8',env:environment}));
}

function fixture(message, initial = false) {
  const dir = mkdtempSync(join(tmpdir(),'multicad-release-contract-'));
  git(dir,['init','-b','main']);
  git(dir,['config','user.name','Release test']);git(dir,['config','user.email','test@example.invalid']);
  writeFileSync(join(dir,'Cargo.toml'),'[workspace.package]\nversion = "0.1.0-dev"\n');
  // A clone has no installed dependencies: preview must bypass its config.
  writeFileSync(join(dir,'release.config.cjs'),"require('multicad-deliberately-uninstalled-fixture');\n");
  git(dir,['add','Cargo.toml','release.config.cjs']);git(dir,['commit','-m','chore: initial fixture']);
  if(!initial)git(dir,['tag','v1.2.3']);
  git(dir,['commit','--allow-empty','-m',message]);
  return dir;
}

for(const [message,expected,initial] of [
  ['fix: repair a reference','1.2.4',false],
  ['perf: improve picking','1.2.4',false],
  ['feat: add a sketch','1.3.0',false],
  ['feat!: replace project format','2.0.0',false],
  ['ci: configure runners',null,false],
  ['feat: first native feature','0.1.0',true],
  ['fix: first reference repair','0.0.1',true],
])test(`actual release-engine preview: ${message}`,async()=>{
  const dir=fixture(message,initial);
  try{
    const before=git(dir,['show-ref']);
    const result=preview(dir);
    assert.equal(result.version,expected);
    assert.equal(git(dir,['show-ref']),before,'preview must not mutate source refs');
    assert.equal(git(dir,['status','--porcelain']),'');
    assert.equal(git(dir,['remote','-v']),'');
    if(result.version)assert.match(result.notes,/sketch|reference|picking|format|feature/i);
  }finally{rmSync(dir,{recursive:true,force:true});}
});

test('an already tagged commit is a no-op on retry',async()=>{
  const dir=fixture('fix: released change');
  try{
    git(dir,['tag','v1.2.4']);const refs=git(dir,['show-ref']);
    assert.equal(preview(dir).version,null);
    assert.equal(git(dir,['show-ref']),refs);
  }finally{rmSync(dir,{recursive:true,force:true});}
});

test('tag-only configuration has no npm or release-commit plugin',()=>{
  assert.deepEqual(config.branches,['main']);assert.equal(config.tagFormat,'v${version}');
  const plugins=config.plugins.map(p=>Array.isArray(p)?p[0]:p);
  assert(!plugins.some(p=>/semantic-release[/\\](?:git|npm)[/\\]/u.test(p)));
  const githubPlugin=require.resolve('@semantic-release/github');
  assert(plugins.includes(githubPlugin));
  const provider=config.plugins.find(p=>Array.isArray(p)&&p[0]===githubPlugin)[1];
  assert.equal(provider.successComment,false);assert.equal(provider.failComment,false);
});

test('publication is disabled unless every trusted activation condition is present',()=>{
  const env={GITHUB_ACTIONS:'true',GITHUB_EVENT_NAME:'push',GITHUB_REF:'refs/heads/main',GITHUB_REPOSITORY:'example/application',GITHUB_SHA:'a'.repeat(40),MULTICAD_EXPECTED_RELEASE_VERSION:'1.2.3',MULTICAD_RELEASES_ENABLED:'true',GITHUB_TOKEN:'test-only'};
  assert.doesNotThrow(()=>requirePublicationEnvironment(env));
  for(const key of Object.keys(env)){const copy={...env};delete copy[key];assert.throws(()=>requirePublicationEnvironment(copy),undefined,key);}
  assert.throws(()=>requirePublicationEnvironment({...env,GITHUB_EVENT_NAME:'pull_request'}));
  assert.throws(()=>requirePublicationEnvironment({...env,GITHUB_REF:'refs/tags/v1.2.3'}));
  assert.throws(()=>githubRepositoryUrl({...env,GITHUB_SERVER_URL:'https://forgejo.example.invalid'}));
});

test('version and repository inputs cannot escape their intended scope',()=>{
  for(const value of ['v1.2.3','1.2','1.2.3\nrelease=true','../1.2.3'])assert.throws(()=>requireVersion(value));
  for(const repository of ['owner/repo/extra','../repo','https://example.invalid/repo'])assert.throws(()=>githubRepositoryUrl({GITHUB_REPOSITORY:repository}));
  assert.equal(githubRepositoryUrl({GITHUB_REPOSITORY:'owner/repo'}),'https://github.com/owner/repo.git');
});

test('native regression suites are read-only and stay runnable by hand',()=>{
  const regression=readFileSync(new URL('../.github/workflows/native-regression.yml',import.meta.url),'utf8');
  const release=readFileSync(new URL('../.github/workflows/release.yml',import.meta.url),'utf8');
  const owned=readFileSync(new URL('../.github/workflows/owned-kernel.yml',import.meta.url),'utf8');
  assert.match(regression,/^permissions:\n {2}contents: read$/m);
  assert.doesNotMatch(regression,/secrets\.|: *write|npm publish|semantic-release|git push|git tag/);
  assert.match(regression,/^ {2}workflow_dispatch:$/m);
  for(const os of ['ubuntu-24\\.04','macos-15','windows-2025'])assert.match(regression,new RegExp(`- os: ${os}\\n`));
  for(const suite of ['--bin multicad files::','--test release_hardening','--test project_file','--test drawing_exports','--test control_desk','cargo nextest run --workspace --locked --no-fail-fast'])assert(regression.includes(suite),suite);
  assert.doesNotMatch(regression,/cargo test /);
  for(const os of ['ubuntu-24\\.04','macos-15','windows-2025'])for(const group of ['kernel','files','interface'])assert.match(regression,new RegExp(`- os: ${os}\\n.*\\n {12}group: ${group}\\n`),`${os} ${group}`);
  // The partitioned complete-workspace archive runs every test these named
  // suites run, so publication no longer waits for a second compile of them.
  assert.doesNotMatch(release,/native-regression\.yml/);
  const packaging=readFileSync(new URL('./package-release.mjs',import.meta.url),'utf8');
  assert.match(packaging,/const testArgs=packageTestArguments\(process\.env\.MULTICAD_PACKAGE_TESTS\?\?'full'\)/);
  assert.match(packaging,/if\(testArgs\)\{\n    run\('cargo',testArgs,\{env,stdio:'inherit'\}\);\n    run\('cargo',workspaceDocTestArguments\(\),\{env,stdio:'inherit'\}\);\n  \}/);
  assert.deepEqual(workspaceDocTestArguments(),['test','--workspace','--doc','--profile','ci','--locked']);
  // Every platform package job runs the packaged example/MCP workflow and the
  // statically linked OpenCascade code audit on the archived executable.
  assert.match(packaging,/rejectOCCTLibraries\(runtimeLibraries\);\n    checkStaticCode\(readFileSync\(layout\.executable\)\);/);
  assert.match(packaging,/await mcpSmoke\(layout\.executable,/);
  const install=`- name: Install pinned cargo-nextest\n        run: cargo install cargo-nextest --version ${nextestVersion} --locked\n`;
  assert.equal(regression.split(install).length-1,1);
  // CI downloads the same pinned nextest release instead of compiling it.
  const prebuilt=`- name: Install pinned cargo-nextest\n        uses: taiki-e/install-action@cfdb7ed4ccc7c7dce37759eae89dc5674b469b39 # v2.90.3\n        with:\n          tool: nextest@${nextestVersion}\n`;
  assert.equal(release.split(prebuilt).length-1,2);
  assert.doesNotMatch(release,/cargo install cargo-nextest/);
  const cargo=readFileSync(new URL('../Cargo.toml',import.meta.url),'utf8');
  const profile=cargo.split('[profile.ci]')[1]?.split(/\n\[/)[0];
  assert(profile,'An optimized CI test profile must exist');
  assert.match(profile,/inherits = "dev"/);
  assert.match(profile,/opt-level = 2/);
  assert.match(profile,/debug-assertions = true/);
  assert.match(profile,/overflow-checks = true/);
  assert.match(owned,/^permissions:\n {2}contents: read$/m);
  for(const os of ['ubuntu-24.04','macos-15','windows-2025'])assert(owned.includes(os));
  assert(owned.includes('cargo test -p multicad-precision -p multicad-rust-kernel --profile ci --locked'));
  assert(owned.includes('cargo clippy -p multicad-precision -p multicad-rust-kernel --all-targets --profile ci --locked -- -D warnings'));
  for(const workflow of [owned,regression,release])assert.doesNotMatch(workflow,/foundation-deps|spikes\/native-foundation|native-deps|features.*occt|MULTICAD_OCCT/);
  assert.match(release,/uses: \.\/\.github\/workflows\/owned-kernel\.yml/);
  // Dependency caches only restore compiled dependencies; every job still
  // checks out and builds the exact commit.
  const cache='uses: Swatinem/rust-cache@d1527c8f67ef9af357290475d4af4b0fab8f12e0 # v2.9.3';
  assert.equal(release.split(cache).length-1,3);
  assert.equal(owned.split(cache).length-1,1);
});

test('active build instructions do not direct developers to the retired provider',()=>{
  for(const file of ['README.md','docs/BUILDING.md','docs/RELEASING.md']){
    const text=readFileSync(new URL('../'+file,import.meta.url),'utf8');
    assert.doesNotMatch(text,/foundation-deps|native-dependencies\.json|target\/native|spikes\/native-foundation|OCCT is statically linked|with exact Open CASCADE/);
    assert.match(text,/Rust/);
  }
});

test('every platform compiles the workspace tests once and its archive runs all gate publication',()=>{
  const release=readFileSync(new URL('../.github/workflows/release.yml',import.meta.url),'utf8');
  const base=['nextest','run','--workspace','--cargo-profile','ci','--locked','--no-fail-fast'];
  assert.deepEqual(workspaceTestArguments(),base);
  assert.deepEqual(workspaceTestArguments('generators'),[...base,'-E',`test(${exampleGeneratorFilter})`]);
  assert.deepEqual(workspaceTestArguments('remaining'),[...base,'-E',`not test(${exampleGeneratorFilter})`]);
  assert.throws(()=>workspaceTestArguments('typo'));
  assert.deepEqual(workspaceListArguments(),['nextest','list','--workspace','--cargo-profile','ci','--locked','--message-format','json']);
  assert.deepEqual(partitionArchiveArguments(['a::x','b::y'],'p.tar.zst'),['nextest','archive','--workspace','--cargo-profile','ci','--locked','-E','binary_id(=a::x) | binary_id(=b::y)','--archive-file','p.tar.zst']);
  assert.throws(()=>partitionArchiveArguments([],'p.tar.zst'));
  assert.deepEqual(archiveRunArguments('p.tar.zst'),['nextest','run','--archive-file','p.tar.zst','--workspace-remap','.','--extract-to','.','--extract-overwrite','--no-fail-fast']);
  assert(Number.isInteger(testPartitions)&&testPartitions>1);
  for(let k=1;k<=testPartitions;k++)assert.equal(partitionArchiveName(k),`nextest-${k}.tar.zst`);
  for(const bad of [0,testPartitions+1,1.5])assert.throws(()=>partitionArchiveName(bad));
  // Every listed binary lands in exactly one partition, whatever the weights
  // say; weights only balance runtime, and the split is deterministic.
  const suites={};
  for(let i=0;i<37;i++)suites[`crate-${i%5}::bin_${i}`]={testcases:Object.fromEntries(Array.from({length:i%7},(_,t)=>[`t${t}`,{}]))};
  const ids=Object.keys(suites).sort();
  for(const weights of [{},{'crate-0::bin_0':10000},{'crate-1::bin_1':0.5,missing:99}]){
    const bins=partitionBinaries(suites,weights);
    assert.equal(bins.length,testPartitions);
    const assigned=bins.flatMap(b=>b.binaries);
    assert.deepEqual([...assigned].sort(),ids);
    assert.equal(new Set(assigned).size,ids.length);
    assert.deepEqual(partitionBinaries(suites,weights),bins);
  }
  const heavy=partitionBinaries(suites,{'crate-0::bin_0':10000});
  assert(heavy.some(b=>b.binaries.length===1&&b.binaries[0]==='crate-0::bin_0'),'the heaviest binary gets its own runner');
  // Local packaging tests first; the CI package job leaves that to the
  // archive runs of the same commit.
  assert.deepEqual(packageTestArguments('full'),base);
  assert.equal(packageTestArguments('gated'),null);
  assert.throws(()=>packageTestArguments('remaining'));
  assert.equal(release.match(/MULTICAD_PACKAGE_TESTS/g)?.length,1);
  assert.match(release,/ {6}MULTICAD_PACKAGE_TESTS: gated\n/);
  const job=name=>release.split(`\n  ${name}:\n`)[1]?.split(/\n  [a-z_]+:\n/)[0]??'';
  const archive=job('test_archive'),test=job('test'),pack=job('package');
  for(const os of ['ubuntu-24.04','macos-15','windows-2025'])for(const [name,text] of [['test_archive',archive],['test',test],['package',pack]])assert(text.includes(`os: ${os}`),`${name} ${os}`);
  assert(archive.includes('run: node tools/native-test-shard.mjs archive nextest'));
  assert(archive.includes('run: cargo test --workspace --doc --profile ci --locked'));
  assert(archive.includes('run: node tools/release-budgets.mjs --profile ci --repeat 1'));
  for(let k=1;k<=testPartitions;k++){
    assert(archive.includes(`name: nextest-\${{ matrix.target }}-${k}-\${{ github.sha }}\n          path: nextest/nextest-${k}.tar.zst\n`),`upload ${k}`);
  }
  assert.equal(archive.match(/path: nextest\/nextest-\d+\.tar\.zst/g)?.length,testPartitions);
  assert.match(test,/needs: test_archive\n/);
  assert(test.includes(`partition: [${Array.from({length:testPartitions},(_,i)=>i+1).join(', ')}]`));
  assert(test.includes('name: nextest-${{ matrix.target }}-${{ matrix.partition }}-${{ github.sha }}'));
  assert(test.includes('run: node tools/native-test-shard.mjs shard "nextest-${{ matrix.partition }}.tar.zst"'));
  assert.doesNotMatch(pack,/nextest/);
  assert.match(release,/needs: \[metadata, quality, package, test, owned-kernel\]/);
  const weights=JSON.parse(readFileSync(new URL('./test-weights.json',import.meta.url),'utf8'));
  assert.match(weights.source,/never skips a test/);
  for(const [id,seconds] of Object.entries(weights.seconds))assert(/^[a-z0-9-]+(::\S+)?$/.test(id)&&Number.isFinite(seconds)&&seconds>=0,id);
});

test('hosted release notes stay below GitHub\'s body limit and link the full change list',async()=>{
  const capped=require('./semantic-release-capped-notes.cjs');
  assert(capped.NOTES_LIMIT<capped.GITHUB_RELEASE_BODY_LIMIT);
  const notesPlugin=config.plugins.find(p=>Array.isArray(p)&&String(p[0]).endsWith('semantic-release-capped-notes.cjs'));
  assert(notesPlugin,'release notes must go through the capped generator');
  assert.deepEqual(notesPlugin[1],{preset:'conventionalcommits'});
  const context={env:{GITHUB_REPOSITORY:'owner/repo'},lastRelease:{gitTag:'v0.1.0'},nextRelease:{gitTag:'v0.2.0'}};
  const short='## 0.2.0\n\n* a change\n';
  assert.equal(capped.capNotes(short,context),short);
  const lines=Array.from({length:6000},(_,i)=>`* change ${i} ${'x'.repeat(40)}`);
  const long=`## 0.2.0\n\n${lines.join('\n')}\n`;
  const result=capped.capNotes(long,context);
  assert(result.length<=capped.NOTES_LIMIT,`${result.length}`);
  assert(result.endsWith('Full list of changes: https://github.com/owner/repo/compare/v0.1.0...v0.2.0\n'));
  const body=result.split('\n\n---\n\n')[0];
  assert(long.startsWith(body),'truncation keeps an unmodified prefix');
  assert(lines.includes(body.split('\n').at(-1)),'truncation keeps whole lines');
  assert.equal(capped.capNotes(long,{},1000).length<=1000,true);
});
