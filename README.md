# Splines
Spline interpolation in Julia


### TODO
# dev
- [x] piecewise interpolation: reduce allocation
- [x] b-splines
  - [x] allow for `closed=true` like piecewise linear
- [x] bezier curve
  - [ ] rational bezier curve
  - [ ] spline by joining multiple d-degree bezier curves
- [ ] algebraic topology to get curve topology from data
- [x] cleanup ploting code
- [x] spline in RGB space
- [ ] parametrize curves somehow
  - [ ] bspline points are not sampled uniformely!
- [ ] add methods to evalute curves at a single parmtere value

# speed
- [ ] why is PWL so slow to optimize first round?
- [X] improve BSpline in-place versino

# tests
- [x] test piecewise linear on 3D and > 3D 

# docs
- [ ] write tests
- [ ] write docs