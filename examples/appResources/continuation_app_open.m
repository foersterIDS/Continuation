clear; close all; clc;

if(~isdeployed)
  cd(fileparts(which('Continuation.mlapp')));
end

open Continuation.mlapp;
