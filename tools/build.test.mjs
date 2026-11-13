import assert from 'node:assert/strict';
import { mkdtempSync, mkdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';
import { assembleDevelopmentRuntime, bootstrap, developmentExecutable, root } from './build.mjs';

test('Rust bootstrap uses only the native toolchain and locked Cargo dependencies',()=>{
  for(const platform of ['darwin','linux','win32']){
    const probes=[],commands=[];
    bootstrap({platform,nodeVersion:'22.14.0',captureCommand:(...call)=>probes.push(call),runCommand:(...call)=>commands.push(call)});
    assert.deepEqual(probes,[['rustup',['--version']],['tar',['--version']],...(platform==='darwin'?[['xcrun',['--find','clang++']]]:[])]);
    assert.deepEqual(commands,[['rustup',['show','active-toolchain']],['cargo',['fetch','--locked']]]);
  }
});

test('bootstrap reports missing Rust and Apple SDK tools before building',()=>{
  for(const [missing,message] of [['rustup',/Install Rust/],['xcrun',/Install Apple Command Line Tools/]]){
    assert.throws(()=>bootstrap({platform:'darwin',nodeVersion:'24.0.0',captureCommand:command=>{if(command===missing)throw new Error('missing');},runCommand:()=>assert.fail('build must not start')}),message);
  }
});

test('failed locked dependency fetch propagates without another bootstrap path',()=>{
  const failure=new Error('offline dependency unavailable'),commands=[];
  assert.throws(()=>bootstrap({platform:'linux',nodeVersion:'24.0.0',captureCommand:()=>{},runCommand:(command,args)=>{commands.push([command,args]);if(command==='cargo')throw failure;}}),error=>error===failure);
  assert.equal(commands.length,2);
});

test('unsupported Node versions fail before probing or fetching',()=>{
  for(const nodeVersion of ['20.19.0','22.13.1'])assert.throws(()=>bootstrap({nodeVersion,captureCommand:()=>assert.fail('unexpected probe'),runCommand:()=>assert.fail('unexpected command')}),/22.14 or newer/);
});

test('development bundles select Rust-only runtime assembly',()=>{
  const calls=[];
  assembleDevelopmentRuntime('/stage','/binary','0.1.0-test',(...call)=>calls.push(call));
  assert.deepEqual(calls,[['/stage','/binary','0.1.0-test',{rustOnly:true}]]);
});

test('development packaging follows the actual Cargo target override',()=>{
  const directory=mkdtempSync(join(tmpdir(),'multicad-cargo-output-'));
  try{
    const expected=join(directory,'debug',process.platform==='win32'?'multicad.exe':'multicad');
    mkdirSync(join(directory,'debug'));
    writeFileSync(expected,'fresh selected output');
    const actual=developmentExecutable({...process.env,CARGO_TARGET_DIR:directory});
    assert.equal(actual,expected);
    assert.equal(readFileSync(actual,'utf8'),'fresh selected output');
  }finally{rmSync(directory,{recursive:true,force:true});}
});

test('normal Cargo and developer packaging agree on the configured output',()=>{
  const environment={...process.env};delete environment.CARGO_TARGET_DIR;
  assert.equal(developmentExecutable(environment),join(root,'target/cargo/debug',process.platform==='win32'?'multicad.exe':'multicad'));
});
