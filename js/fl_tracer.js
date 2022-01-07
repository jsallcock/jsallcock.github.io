// EQVIS
// J. S. Allcock 2021
// 
// At its core a JavaScript port of Nick Walkden's pyFieldlineTracer. 
// I am *NOT* experienced with JavaScript so this could cause big problemos (memory leaks etc.)

// import test data: normalised magnetic field vector (R, Z, phi) on cylindrical coordinate grid (R, Z) (toroidal symmetry assumed) 
import data from './45245.js';
var IDX_TIME = 0;
var ONE_SIXTH = 1 / 6;
var DEFAULT_MAXSTEPS = 3000;
var DEFAULT_DS = 1e-2;
var DEFAULT_RSTART = 1.3;
var DEFAULT_NLINES = 5;
var DEFAULT_TOR_LIM = 2 * 3.1415;
var DEFAULT_LINEWIDTH = 1.5;
var DEFAULT_SHOWAXIS = false;
var DEFAULT_SHOWXPOINTS = false;
var PI = 3.14159265359;

//-------------------------------------------------------------------------------
// Tools
//-------------------------------------------------------------------------------
class Blerp {

  constructor (x, y, z) {
    // Bilinear interpolation, benchmarked against scipy.interpolate.interp2d.
    // x is column coords and y is row coords (matching scipy.interpolate.interp2d convention).
    // x and y must be monotonic arrays.
    // z is an array of arrays, indexed like z[idx_row][idx_column]
    this.x = x;
    this.y = y;
    this.z = z;
  }

  call(x, y) {
    // Interpolate z at given (x, y). 
    // Returns NaN if x, y is outside original xy grid.
    let x2_idx = this.x.findIndex(n => n > x);
    let x1_idx = x2_idx - 1;
    let y2_idx = this.y.findIndex(n => n > y);
    let y1_idx = y2_idx - 1;
    if (x1_idx < 0 || y1_idx < 0 || x2_idx == -1 || y2_idx == -1) { 
      return NaN 
    };
    let x2 = this.x[x2_idx];
    let x1 = this.x[x1_idx];
    let y2 = this.y[y2_idx];
    let y1 = this.y[y1_idx];
    let q11 = this.z[y1_idx][x1_idx];
    let q12 = this.z[y2_idx][x1_idx];
    let q21 = this.z[y1_idx][x2_idx];
    let q22 = this.z[y2_idx][x2_idx];
    let norm = 1 / ((x2 - x1) * (y2 - y1));
    let x2mx = x2 - x;
    let xmx1 = x - x1;
    let y2my = y2 - y;
    let ymy1 = y - y1;
    let w11 = x2mx * y2my * norm;
    let w12 = x2mx * ymy1 * norm;
    let w21 = xmx1 * y2my * norm;
    let w22 = xmx1 * ymy1 * norm;
    return w11 * q11 + w12 * q12 + w21 * q21 + w22 * q22
  }
}

function cyl2cart(r, z, phi) { return [r * Math.cos(phi), r * Math.sin(phi), z] }

var r_wall = [1.56442 , 1.73298 , 1.34848 , 1.0882  , 0.902253, 0.903669,
  0.533866, 0.538011, 0.332797, 0.332797, 0.334796, 0.303115,
  0.305114, 0.269136, 0.271135, 0.260841, 0.260841, 0.271135,
  0.269136, 0.305114, 0.303115, 0.334796, 0.332797, 0.332797,
  0.538598, 0.534469, 0.90563 , 0.904219, 1.0882  , 1.34848 ,
  1.73018 , 1.56442 , 1.37999 , 1.37989 , 1.19622 , 1.19632 ,
  1.05537 , 1.05528 , 0.947502, 0.905686, 0.899143, 0.883388,
  0.867681, 0.851322, 0.833482, 0.826063, 0.822678, 0.821023,
  0.820691, 0.822887, 0.827573, 0.839195, 0.855244, 0.877567,
  0.899473, 1.18568 , 1.279   , 1.296   , 1.521   , 1.521   ,
  1.8     , 1.8     , 1.521   , 1.521   , 1.296   , 1.279   ,
  1.18568 , 0.899473, 0.877567, 0.855244, 0.839195, 0.827573,
  0.822887, 0.820691, 0.821023, 0.822678, 0.826063, 0.833482,
  0.851322, 0.867681, 0.883388, 0.899143, 0.905686, 0.947502,
  1.05528 , 1.05537 , 1.19632 , 1.19622 , 1.37989 , 1.37999 ,
  1.56442 ];
var z_wall = [ 1.56424 ,  1.67902 ,  2.06041 ,  2.05946 ,  1.87565 ,  1.87424 ,
  1.50286 ,  1.49874 ,  1.29709 ,  1.094   ,  1.094   ,  0.8475  ,
  0.8475  ,  0.565   ,  0.565   ,  0.495258, -0.507258, -0.577   ,
 -0.577   , -0.8595  , -0.8595  , -1.106   , -1.106   , -1.30909 ,
 -1.5099  , -1.51403 , -1.88406 , -1.88547 , -2.06614 , -2.06519 ,
 -1.68099 , -1.56884 , -1.57688 , -1.57673 , -1.58475 , -1.5849  ,
 -1.59105 , -1.59091 , -1.59561 , -1.59556 , -1.59478 , -1.59026 ,
 -1.58087 , -1.56767 , -1.54624 , -1.52875 , -1.51517 , -1.49624 ,
 -1.47724 , -1.44582 , -1.41923 , -1.38728 , -1.35284 , -1.3221  ,
 -1.30018 , -1.0138  , -0.8423  , -0.8202  , -0.8202  , -0.25    ,
 -0.25    ,  0.25    ,  0.25    ,  0.8156  ,  0.8156  ,  0.8377  ,
  1.0092  ,  1.29558 ,  1.3175  ,  1.34824 ,  1.38268 ,  1.41463 ,
  1.44122 ,  1.47264 ,  1.49164 ,  1.51057 ,  1.52415 ,  1.54164 ,
  1.56307 ,  1.57627 ,  1.58566 ,  1.59018 ,  1.59096 ,  1.59101 ,
  1.58631 ,  1.58645 ,  1.5803  ,  1.58015 ,  1.57213 ,  1.57228 ,
  1.56424 ];

//-------------------------------------------------------------------------------
// Two-line intersection check tools
//-------------------------------------------------------------------------------
// Source: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
class Point {

  constructor(r, z) {
    this.r = r;
    this.z = z;
  }
}

// Given three collinear points p1, p2, p3, the function checks if
// point p2 lies on line segment p1 --> p3
function onSegment(p1, p2, p3){
  if (p2.r <= Math.max(p1.r, p3.r) && p2.r >= Math.min(p1.r, p3.r) &&
      p2.z <= Math.max(p1.z, p3.z) && p2.z >= Math.min(p1.z, p3.z)) return true;  
  return false;
}

// Find orientation of ordered triplet (p1, p2, p3).
// function returns following values:
// 0 --> p1, p2 and p3 are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
function orientation(p1, p2, p3) {
  let val = (p2.z - p1.z) * (p3.r - p2.r) - (p2.r - p1.r) * (p3.z - p2.z);
  if (val == 0) return 0;
  return (val > 0)? 1: 2;
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
function testIntersect(p1, q1, p2, q2) {
    // Find the four orientations needed for general and special cases
    let o1 = orientation(p1, q1, p2);
    let o2 = orientation(p1, q1, q2);
    let o3 = orientation(p2, q2, p1);
    let o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
    return false; // Doesn't fall in any of the above cases
}

//-------------------------------------------------------------------------------
// Field line tracer
//-------------------------------------------------------------------------------
class Tracer {
  constructor (canvas_width, canvas_height) {
    let r = data.r;
    let z = data.z;
    this.canvas_width = canvas_width;
    this.canvas_height = canvas_height;
    this.br = new Blerp(r, z, data.Br);
    this.bz = new Blerp(r, z, data.Bz);
    this.bphi = new Blerp(r, z, data.Bphi);
    this.ds = DEFAULT_DS;

    this.rstart = DEFAULT_RSTART;
    this.nlines = DEFAULT_NLINES;
    this.maxSteps = DEFAULT_MAXSTEPS;
    this.showAxis = DEFAULT_SHOWAXIS;
    this.showXpoints = DEFAULT_SHOWXPOINTS;

    // this is terrible 
    this.rstartLastUsed = null;
    this.nlinesLastUsed = null;
    this.maxStepsLastUsed = null;
    this.showAxisLastUsed = null;
    this.showXpointsLastUsed = null;

    this.material = new THREE.LineBasicMaterial( {
      color: 0xffffff,
      linewidth: DEFAULT_LINEWIDTH,
      linejoin: 'round'
    } );

    // initialise
    this.scene = new THREE.Scene();
    this.lines = [];
    this.init();
    this.update();
    this.update_axis();
    this.update_xpoints();
  }

  init () {
    // initialise buffers for the field lines
    // this must be called whenever nlines changes
    // following: https://threejs.org/docs/#manual/en/introduction/How-to-update-things
    if (this.nlinesLastUsed != null) {
      for (let i = 0; i < this.nlinesLastUsed; i++) { this.remove_from_scene( 'fieldline' + i.toString() ) }
    }
    this.lines = []
    for (let i = 0; i < this.nlines; i++) {
      const geometry = new THREE.BufferGeometry();
      // const positions = new Float64Array( N_STEPS_MAX * 3 ); // 3 vertices per point
      const positions = new Float32Array( this.maxSteps * 3 ); // 3 vertices per point
      geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
      const line = new THREE.Line( geometry,  this.material );
      line.name = 'fieldline' + i.toString();
      this.lines.push( line );
      this.scene.add( line );   
    }
    this.nlinesLastUsed = this.nlines;
    this.maxStepsLastUsed = this.maxSteps;
  }

  update () {
    // re-trace field lines and update scene
    // (actually just traces one and adds the appropriate toroidal angle!)
    let [r, z, phi] = this.trace( this.rstart, 0, 0, );
    for (let j = 0; j < this.nlines; j++) {
      let delta_phi = j * 2 * PI / this.nlines;
      let positions = this.lines[ j ].geometry.attributes.position.array;
      let idx = 0;
      for (let i = 0; i < r.length; i++) {
        if (!isNaN(r[i])) {
          let [x_i, y_i, z_i] = cyl2cart(r[i], z[i], phi[i] + delta_phi )
          positions[ idx ++ ] = x_i;
          positions[ idx ++ ] = y_i;
          positions[ idx ++ ] = z_i;
        }
      }
      this.lines[ j ].geometry.setDrawRange( 0, r.length ); // only need to draw the points that are there
      this.lines[ j ].geometry.attributes.position.needsUpdate = true;
      this.lines[ j ].geometry.computeBoundingBox();
      this.lines[ j ].geometry.computeBoundingSphere();
    }
    this.rstartLastUsed = this.rstart;
  }

  update_axis () {
    if (this.showAxis && this.scene.getObjectByName('axis') == undefined) {
      this.addCircle(data['r_axis'], data['z_axis'], 'Crimson', 'axis' );
    }

    if (!this.showAxis && this.scene.getObjectByName('axis') != undefined) {
      this.remove_from_scene('axis');
    }
    this.showAxisLastUsed = this.showAxis;
  }

  update_xpoints () {
    if (this.showXpoints && this.scene.getObjectByName('xpoints') == undefined) {
      this.addCircle(data['lower_xpoint_r'], data['lower_xpoint_z'], 'Chartreuse', 'xpoint_lower' );
      this.addCircle(data['upper_xpoint_r'], data['upper_xpoint_z'], 'Chartreuse', 'xpoint_upper' );
    }
    if (!this.showXpoints && this.scene.getObjectByName('xpoint_upper') != undefined) {
      this.remove_from_scene('xpoint_lower');
      this.remove_from_scene('xpoint_upper');
    }
    this.showXpointsLastUsed = this.showXpoints;
  }

  trace ( rstart, zstart, phistart ) {
    // trace single fieldline, returning array of r, z, phi coords
    let [r_fwd, z_fwd, phi_fwd] = this.trace1way(rstart, zstart, phistart, this.ds);
    let [r_bck, z_bck, phi_bck] = this.trace1way(rstart, zstart, phistart, -this.ds);
    let r = r_bck.reverse().concat(r_fwd);
    let z = z_bck.reverse().concat(z_fwd);
    let phi = phi_bck.reverse().concat(phi_fwd);
    return [r, z, phi]
  }

  trace1way ( rstart, zstart, phistart, ds ) {
    // trace single fieldline one-way, returning array of r, z, phi coords
    // TODO check toroidal angle limit
    let r = [rstart, ];
    let z = [zstart, ];
    let phi = [phistart, ];
    let i = 0;
    while (i < this.maxSteps / 2) {  
      // factor of 1 / 2 since we trace both ways. maxSteps refers to composite field line. 
      // This is important because I use maxSteps to create buffers for the field lines.
      const [r_i, z_i, phi_i] = [r[i], z[i], phi[i]];
      const [dr, dz, dphi] = this.take_step(r_i, z_i, ds);
      const [r_ip, z_ip, phi_ip] = [r_i + dr, z_i + dz, phi_i + dphi];
      if (this.testIntersectWall(r_i, z_i, r_ip, z_ip)) {
        break;
      }
      r.push(r_ip);
      z.push(z_ip);
      phi.push(phi_ip);
      i++;
    }
    return [r, z, phi]
  }

  take_step ( r, z, ds, ) {
    // Runge-Kutta method (RK4) to take a single step 
    // (r, z) are coordinates of starting point (m), ds is step size (m)
    // NOT YET BENCHMARKED
    let dr1 = ds * this.br.call(r, z);
    let dz1 = ds * this.bz.call(r, z);
    let dphi1 = ds * this.bphi.call(r, z) / r;
    let dr1_half = 0.5 * dr1;
    let dz1_half = 0.5 * dz1;
    let dr2 = ds * this.br.call(r + dr1_half, z + dz1_half);
    let dz2 = ds * this.bz.call(r + dr1_half, z + dz1_half);
    let dphi2 = ds * this.bphi.call(r + dr1_half, z + dz1_half) / r;
    let dr2_half = 0.5 * dr2;
    let dz2_half = 0.5 * dz2;
    let dr3 = ds * this.br.call(r + dr2_half, z + dz2_half);
    let dz3 = ds * this.bz.call(r + dr2_half, z + dz2_half);
    let dphi3 = ds * this.bphi.call(r + dr2_half, z + dz2_half) / r;
    let dr4 = ds * this.br.call(r + dr3, z + dz3);
    let dz4 = ds * this.bz.call(r + dr3, z + dz3);
    let dphi4 = ds * this.bphi.call(r + dr3, z + dz3) / r;
    let dr = ONE_SIXTH * (dr1 + 2 * dr2 + 2 * dr3 + dr4); 
    let dz = ONE_SIXTH * (dz1 + 2 * dz2 + 2 * dz3 + dz4);
    let dphi = ONE_SIXTH * (dphi1 + 2 * dphi2 + 2 * dphi3 + dphi4);
    return [dr, dz, dphi]
  }

  testIntersectWall ( r1, z1, r2, z2 ) {
    // does the line segment (r1, z1) --> (r2, z2) intersect with the MAST-U wall?
    // this currently leaves a small gap through which a field line could pass -- TODO: fix
    for (let i = 0; i < r_wall.length - 1; i++) {
      let p1_fl = new Point(r1, z1);
      let p2_fl = new Point(r2, z2);
      let p1_w = new Point(r_wall[i], z_wall[i]);
      let p2_w = new Point(r_wall[i+1], z_wall[i+1]);
      if ( testIntersect(p1_fl, p2_fl, p1_w, p2_w) ) return true;
    }
    return false; 
  }

  addCircle (r, z_centre, color='Crimson', name='') {
    // Circle parallel to x-y plane
    const segmentCount = 100;
    let points = []
    for (var i = 0; i <= segmentCount; i++) {
      let theta = (i / segmentCount) * Math.PI * 2;
      points.push( new THREE.Vector3( Math.cos(theta) * r, Math.sin(theta) * r, z_centre ));            
    }
    // const material = new THREE.LineBasicMaterial({ color: color, linewidth: DEFAULT_LINEWIDTH});
    // const geometry = new THREE.BufferGeometry().setFromPoints( points );
    // this.scene.add( new THREE.Line( geometry, material ) );
    const line = new MeshLine();
    line.setPoints(points);
    const material = new MeshLineMaterial(
      { lineWidth: 0.025,
        resolution: new THREE.Vector2( this.canvas_width, this.canvas_height ),
        color: color,
      }
    );
    const mesh = new THREE.Mesh(line, material);
    mesh.name = name;
    this.scene.add(mesh);
  }
  remove_from_scene ( name ) {
    this.scene.getObjectByName( name ).geometry.dispose();
    this.scene.getObjectByName( name ).material.dispose();
    this.scene.remove( this.scene.getObjectByName( name ) );
  }
}

//-------------------------------------------------------------------------------
// Make animation
//-------------------------------------------------------------------------------
function plot_fl(){
  let cam, scene, renderer, geometry, material, controls, line;
  var canvas = document.getElementById('canvas-holder-fltracer');
  var width = canvas.clientWidth;
  var height = width * 0.35;
  var tracer = new Tracer(width, height);
  renderer = new THREE.WebGLRenderer({ antialias: true } );
  renderer.setSize( width, height );
  canvas.appendChild( renderer.domElement );
  init();
  animate();

  function init() {   
    cam = new THREE.PerspectiveCamera( 60, width / height, 0.01, 500 );  // args: fov, aspect, near plane, far plane
    cam.position.set( 3, 0, 0 );
    cam.up.set(0, 0, 1)
    cam.lookAt( 0, 0, 0 );
    
    const panel = new dat.GUI( { width: 310, autoPlace: false } );
    document.getElementById('dat-gui-holder-fltracer').appendChild(panel.domElement);
    panel.add(tracer, 'rstart', 0.9, 1.5, 0.005);
    panel.add(tracer, 'nlines', 1, 20, 1);
    panel.add(tracer, 'maxSteps', 800, 10000, 200);
    panel.add(tracer, 'showAxis');
    panel.add(tracer, 'showXpoints');
    
    controls = new THREE.OrbitControls(cam, renderer.domElement);
    controls.update();
  }
  function animate() {
    requestAnimationFrame( animate );
    controls.update();

    if (tracer.nlines != tracer.nlinesLastUsed) {
      renderer.renderLists.dispose();
      tracer.init()
      tracer.update()
    }
    if (tracer.rstart != tracer.rstartLastUsed) {
      tracer.update()
    }
    if (tracer.maxSteps != tracer.maxStepsLastUsed) {
      tracer.init()
      tracer.update()
    }
    if (tracer.showAxis != tracer.showAxisLastUsed) {
      tracer.update_axis()
    }
    if (tracer.showXpoints != tracer.showXpointsLastUsed) {
      tracer.update_xpoints()
    }

    renderer.render( tracer.scene, cam );
  }
}
plot_fl();
// console.log(bf_data);