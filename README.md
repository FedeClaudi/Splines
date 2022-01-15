# Splines
Spline interpolation in Julia


### TODO
# dev
- [x] piecewise interpolation: reduce allocation
- [x] b-splines
  - [x] allow for `closed=true` like piecewise linear
- [x] bezier curve
  - [x] rational bezier curve
  - [ ] spline by joining multiple d-degree bezier curves
- [ ] algebraic topology to get curve topology from data
- [x] cleanup ploting code
- [x] spline in RGB space
- [x] parametrize curves somehow
  - [ ] piecewise linear is not parametrized!
- [x] add methods to evalute curves at a single parmtere value
- [ ] compute curve geometry (length, curvature etc...)
- [ ] consider mergin Bezier with RationalBezier and BSpline with rational BSpline (and make rational BSpline)

- [ ] Two dimensional splines
  - [ ] Bezier surface
  - [ ] NURBS

- [ ] Generalized n-dim splines

# speed
- [ ] why is PWL so slow to optimize first round?
- [X] improve BSpline in-place versino

# tests
- [x] test piecewise linear on 3D and > 3D 

# docs
- [ ] write tests
- [ ] write docs
- [ ] release