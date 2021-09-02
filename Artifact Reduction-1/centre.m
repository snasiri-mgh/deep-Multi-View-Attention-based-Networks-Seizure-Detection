function [sc,mm]=centre(s)
[n,N]=size(s);

sc=s-diag(mean(s'))*ones(n,N);

mm=(mean(s'));