import assert from 'node:assert/strict';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import test from 'node:test';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import { archiveName, sourceName, digest, targets, validateDescriptor, validateDiagnostics, verifyPackageSet, verifyPackage, zip } from './package-release.mjs';

const version='1.2.3',commit='a'.repeat(40);
function fixture(){
  const dir=mkdtempSync(join(tmpdir(),'multicad-package-contract-'));
  const stage=mkdtempSync(join(tmpdir(),'multicad-package-stage-'));
  try{
    for(const target of targets){
      const identity={version,tag:`v${version}`,commit,target,kind:'native-development-host',geometryBackend:'rust'};
      const name=archiveName(version,target),sources=sourceName(version,target);
      for(const asset of [name,sources]){
        const stem=asset.slice(0,-4);mkdirSync(join(stage,stem));
        writeFileSync(join(stage,stem,'BUILD-IDENTITY.json'),JSON.stringify(identity));
        if(asset===name){
          const [os,arch]=target.split('-');
          writeFileSync(join(stage,stem,'BUILD-IDENTITY.txt'),`ToskCAD ${version}\nTarget: ${os}/${{x64:'x86_64',arm64:'aarch64'}[arch]}\nSource commit: ${commit}\nCAD kernel: MultiCAD Rust\n`);
        }
        zip(stage,stem,join(dir,asset));
        writeFileSync(join(dir,asset+'.sha256'),`${digest(join(dir,asset))}  ${asset}\n`);
      }
      writeFileSync(join(dir,name+'.json'),JSON.stringify({...identity,archive:name,sha256:digest(join(dir,name)),sources:{archive:sources,sha256:digest(join(dir,sources))},runtimeLibraries:[]}));
    }
    return dir;
  }finally{rmSync(stage,{recursive:true,force:true});}
}

test('aggregate verification requires all three exact version/commit/target identities and is repeatable',()=>{
  const dir=fixture();try{
    const manifest=verifyPackageSet(dir,version,commit);assert.equal(manifest.assets.length,3);assert.equal(manifest.geometryBackend,'rust');
    writeFileSync(join(dir,'release-manifest.json'),JSON.stringify(manifest));
    assert.deepEqual(verifyPackageSet(dir,version,commit),manifest);
  }finally{rmSync(dir,{recursive:true,force:true});}
});
for(const scenario of ['missing target','changed bytes','wrong commit','wrong tag','stale asset','missing sources','changed sources'])test(`package validation rejects ${scenario}`,()=>{
  const dir=fixture();try{
    const file=join(dir,archiveName(version,targets[0]));
    if(scenario==='missing target')rmSync(file);
    if(scenario==='changed bytes')writeFileSync(file,'damaged');
    if(scenario==='wrong commit'){const d=JSON.parse(readFileSync(file+'.json'));d.commit='b'.repeat(40);writeFileSync(file+'.json',JSON.stringify(d));}
    if(scenario==='wrong tag'){const d=JSON.parse(readFileSync(file+'.json'));d.tag='v9.9.9';writeFileSync(file+'.json',JSON.stringify(d));}
    if(scenario==='stale asset')writeFileSync(join(dir,'old-release.zip'),'stale');
    if(scenario==='missing sources')rmSync(join(dir,sourceName(version,targets[0])));
    if(scenario==='changed sources')writeFileSync(join(dir,sourceName(version,targets[0])),'changed');
    assert.throws(()=>verifyPackageSet(dir,version,commit));
  }finally{rmSync(dir,{recursive:true,force:true});}
});

test('descriptor validation rejects an unknown target and invalid SemVer',()=>{
  assert.throws(()=>validateDescriptor({},version,commit,'unknown','hash'));
  assert.throws(()=>archiveName('../bad','linux-x64'));
});

const releaseTarget='macos-arm64';
const releaseIdentity={version,commit,target:releaseTarget,tag:`v${version}`,kind:'native-development-host',geometryBackend:'rust'};
const releaseDiagnostics=`ToskCAD ${version}\nTarget: macos/aarch64\nSource commit: ${commit}\nCAD kernel: MultiCAD Rust\n`;

test('release diagnostics demand exact version, commit, target and Rust-only kernel',()=>{
  const verify=text=>validateDiagnostics(text,version,commit,releaseTarget);
  verify(releaseDiagnostics);
  verify(releaseDiagnostics.replaceAll('\n','\r\n'));
  for(const text of [
    releaseDiagnostics.replace(version,'1.2.30'),
    releaseDiagnostics.replace(commit,'b'.repeat(40)),
    releaseDiagnostics.replace('macos/aarch64','macos/x86_64'),
    releaseDiagnostics.replace('MultiCAD Rust','OCCT 7.9.3 (statically linked)'),
    releaseDiagnostics+'OCCT source: stale\n',
    releaseDiagnostics+'CAD kernel: OCCT 7.9.3 (statically linked)\n',
  ])assert.throws(()=>verify(text),/diagnostics/);
});

test('release archive verification binds both ZIPs, binary diagnostics and the release descriptor',()=>{
  const dir=mkdtempSync(join(tmpdir(),'multicad-rust-release-'));
  try{
    const name=archiveName(version,releaseTarget);
    const source=name.replace(/\.zip$/,'-sources.zip');
    const runtimeIdentity=join(dir,name.slice(0,-4),'BUILD-IDENTITY.json');
    const runtimeDiagnostics=join(dir,name.slice(0,-4),'BUILD-IDENTITY.txt');
    const sourceIdentity=join(dir,source.slice(0,-4),'BUILD-IDENTITY.json');
    for(const stem of [name,source].map(file=>file.slice(0,-4)))mkdirSync(join(dir,stem));
    writeFileSync(runtimeIdentity,JSON.stringify(releaseIdentity));
    writeFileSync(runtimeDiagnostics,releaseDiagnostics);
    writeFileSync(sourceIdentity,JSON.stringify(releaseIdentity));
    const adapter=join(dir,source.slice(0,-4),'crates/multicad-ui/src');
    mkdirSync(adapter,{recursive:true});
    writeFileSync(join(adapter,'native_window.mm'),'#import <AppKit/AppKit.h>\n#import <Foundation/Foundation.h>\n');
    writeFileSync(join(adapter,'native_window.hpp'),'// retained macOS window adapter header\n');
    const archive=asset=>{
      const file=join(dir,asset);rmSync(file,{force:true});zip(dir,asset.slice(0,-4),file);
      const hash=digest(file);writeFileSync(file+'.sha256',`${hash}  ${asset}\n`);return hash;
    };
    const descriptor={...releaseIdentity,archive:name,sha256:archive(name),sources:{archive:source,sha256:archive(source)},runtimeLibraries:['/usr/lib/libSystem.B.dylib']};
    const writeDescriptor=()=>writeFileSync(join(dir,name+'.json'),JSON.stringify(descriptor));
    const verify=()=>verifyPackage(dir,version,commit,releaseTarget);
    writeDescriptor();
    assert.deepEqual(verify(),descriptor);
    assert.doesNotThrow(()=>validateDescriptor(descriptor,version,commit,releaseTarget,descriptor.sha256));

    descriptor.runtimeLibraries.push('/usr/lib/libTKernel.dylib');
    writeDescriptor();assert.throws(verify,/OpenCascade/);
    descriptor.runtimeLibraries.pop();
    descriptor.nativeSourceHash='c'.repeat(64);
    writeDescriptor();assert.throws(verify,/release identity/);
    delete descriptor.nativeSourceHash;

    writeFileSync(sourceIdentity,JSON.stringify({...releaseIdentity,commit:'b'.repeat(40)}));
    descriptor.sources.sha256=archive(source);
    writeDescriptor();assert.throws(verify,/release identity/);
    writeFileSync(sourceIdentity,JSON.stringify(releaseIdentity));
    descriptor.sources.sha256=archive(source);

    // Each fixture is an actual ZIP with updated valid hashes and identities.
    // Removing a retired artifact must restore admission without weakening the
    // manifest or binary identity contract.
    for(const [asset,path] of [
      [source,'artifacts/native-deps/occt-7.9.3.tar.gz'],
      [source,'crates/multicad-native-build/src/lib.rs'],
      [source,'crates/multicad-kernel/src/bridge.cpp'],
      [source,'crates/multicad-kernel/assets/thread-prototypes/m8.brep'],
      [source,'tools/generate-thread-prototypes.mjs'],
      [name,'resources/OCCT-LGPL-EXCEPTION.txt'],
      [name,'resources/native-dependencies.json'],
      [name,'lib/libTKernel.dylib'],
    ]){
      const leaf=join(dir,asset.slice(0,-4),path);
      const created=[];
      for(let parent=dirname(leaf);!existsSync(parent);parent=dirname(parent))created.push(parent);
      mkdirSync(dirname(leaf),{recursive:true});writeFileSync(leaf,'retired fixture');
      if(asset===name)descriptor.sha256=archive(asset);else descriptor.sources.sha256=archive(asset);
      writeDescriptor();assert.throws(verify,/retired geometry asset/);
      rmSync(leaf);for(const directory of created)rmSync(directory,{recursive:true});
      if(asset===name)descriptor.sha256=archive(asset);else descriptor.sources.sha256=archive(asset);
      writeDescriptor();assert.deepEqual(verify(),descriptor);
    }

    writeFileSync(runtimeDiagnostics,releaseDiagnostics.replace('MultiCAD Rust','OCCT'));
    descriptor.sha256=archive(name);
    writeDescriptor();assert.throws(verify,/diagnostics/);
    writeFileSync(runtimeDiagnostics,releaseDiagnostics);
    descriptor.sha256=archive(name);
    writeDescriptor();
    writeFileSync(join(dir,source+'.sha256'),'stale checksum');
    assert.throws(verify,/Checksum mismatch/);
  }finally{rmSync(dir,{recursive:true,force:true});}
});

test('retired packaging CLI modes cannot select a second build provider',()=>{
  for(const action of ['build-rust-candidate','verify-rust-candidate']){
    const result=spawnSync(process.execPath,[fileURLToPath(new URL('./package-release.mjs',import.meta.url)),action,version],{encoding:'utf8'});
    assert.equal(result.status,1);
    assert.match(result.stderr,/Usage: package-release.mjs build VERSION \| verify VERSION COMMIT/);
  }
});
