# extrap
Simple CBS extrapolation &amp; utility script for creating focal point tables.

## Install
Append to your .bashrc, where `PATH_TO_EXTRAP` is something like `/home/mmd01986/dev/`, then,

```export PYTHONPATH=$PYTHONPATH:PATH_TO_EXTRAP```

## How to use
Example using Feller exponential extrapolation (SCF 3 point extrapolation):

```
>>> import extrap
>>> xs = [2,3,4]
>>> es = [-1.0,-1.1,-1.11]
>>> ans = extrap.extrapolate(xs,es,extrap.exp,pts=3)
>>> ans[0][0]
-1.1111111111111112
>>> 
```

Example using X^-3 (correlation 2 point extrapolation)

```
>>> xs = [3,4]
>>> es = [-1.1,-1.11]
>>> ans = extrap.extrapolate(xs,es,extrap.xm3)
>>> ans[0][0]
-1.1172972972972974
```
