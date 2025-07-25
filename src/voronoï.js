
// helpers for triangulation (combination)
Array.prototype.combine = function (k) {const output = []; // ---> array of values (Array.prototype)
  function backtrack (start=0, path=[]) {
    if (path.length === k) {output.push([...path]); return;};
    for (let i = start; i < this.length; i++) {path.push(this[i]); backtrack.bind(this)(i + 1, path); path.pop();};};
  backtrack.bind(this)(); return output;
};

const ndIndex = (n, k) => {let output = []; // ---> array of indices (Function)
  const backtrack = (start=0, path=[]) => {
    if (path.length === k) {output.push([...path]); return;};
    for (let i = start; i < n; i++) {path.push(i); backtrack(i + 1, path); path.pop();};};
  backtrack(); return output;
};

class Point {
  constructor(x, y) {this._x = x; this._y = y;}
  get x() {return this._x;}
  get y() {return this._y;}
  set x(v) {this._x = v;}
  set y(v) {this._y = v;}
  to(point) {return Math.sqrt((this.x - point.x)**2 + (this.y - point.y)**2);}
};

class Mediatrice { // euclidean equation with shape: ax+by=c
  constructor(p, q) {this._a = p.x - q.x; this._b = p.y - q.y; this._c = (p.x - q.x) * (p.x + q.x) / 2 + (p.y - q.y) * (p.y + q.y) / 2;}

  get a() {return this._a;}
  get b() {return this._b;}
  get c() {return this._c;}
  set a(v) {this._a = v;}
  set b(v) {this._b = v;}
  set c(v) {this._c = v;}

  intercept(line) {
    let {a, b, c} = this, {a:a_, b:b_, c:c_} = line;
    let determinant = a * b_ - a_ * b;
    if (Math.abs(determinant) < 1e-6) {return null;}; // numerical stability
    let [x, y] = [(c * b_ - c_ * b) / determinant, (a * c_ - a_ * c) / determinant];
    return new Point(x, y);}

  intersect(center, radius) {// general quadratic case formulas
    if (this.b !== 0) { // non vertical line
      let A = 1 + (this.a / this.b)**2;
      let B = 2 * (-this.c * this.a / this.b**2 - center.x + this.a / this.b * center.y);
      let C = center.x**2 + center.y**2 - (2 * this.c / this.b * center.y) + (this.c / this.b)**2 - radius**2;
      let discriminant = B**2 - 4 * A * C; // discriminant is always strictly positive (two solutions)
      let x1 = (-B + Math.sqrt(discriminant)) / (2 * A); let y1 = (- this.a * x1 + this.c) / this.b;
      let x2 = (-B - Math.sqrt(discriminant)) / (2 * A); let y2 = (- this.a * x2 + this.c) / this.b;
      return [new Point(x1, y1), new Point(x2, y2)];
    } else { // vertical line
      let B = -2 * center.y;
      let C = center.x**2 + center.y**2 - 2 * this.c / this.a * center.x + (this.c / this.a)**2 - radius**2;
      let discriminant = B**2 - 4*C; // discriminant is always strictly positive (two solutions)
      let y1 = (-B + Math.sqrt(discriminant)) / 2; let y2 = (-B - Math.sqrt(discriminant)) / 2;
      return [new Point(this.c / this.a, y1), new Point(this.c / this.a, y2)];};}
};

class Voronoï extends Map {
  constructor(array) {super(); this.setSeeds(array); this._borders = this.setBorders(array); this._mediatrices = this.setMediatrices();}
  setSeeds(array) {
    // array.sort((p,q) => p.y - q.y); array.sort((p,q) => p.x - q.x); ---> For future algorithm improvement
    for (let point of array) {this.set(point, {vertices: [], edges: {}});};}
  setBorders(array) {// set borders equations
    let U = 0, D = 0, L = 0, R = 0;
    for (let {x, y} of array) {L = Math.min(L, x); R = Math.max(R, x); U = Math.max(U, y); D = Math.min(D, y);};
    L--; R++; D--; U++;
    return {L, R, U, D};}
  setMediatrices() {// unique perpendicular bissectors
    let lines = {}, seeds = [...this.keys()];
    for (let i = 0; i < this.size; i++) {for (let j = i + 1; j < this.size; j++) {lines[`${i}${j}`] = new Mediatrice (seeds[i], seeds[j]);};};
    return lines;}
  cells() {// Delaunay triangulation
    let seeds = [...this.keys()], combinations = ndIndex(this.size, 3);
    for (let triangle of combinations) {
      let segments = [`${triangle[0]}${triangle[1]}`, `${triangle[1]}${triangle[2]}`, `${triangle[0]}${triangle[2]}`];
      let mediatrices = segments.map(side => this._mediatrices[side]);
      let center = mediatrices[0].intercept(mediatrices[1]);
      if (!center) {continue;}; // colinear points --> Delaunay triangulation degenerates (flat triangle)
      let radius = seeds[triangle[0]].to(center);
      if (seeds.filter((_, i) => triangle.includes(i) === false).some(p => p.to(center) < radius)) {continue;}; // Delaunay condition not matched (has nearest neighbour)
      let coCirculars = seeds.filter((p, i) => triangle.includes(i) === false && p.to(center) === radius);
      triangulation: for (let k = 0; k < 3; k++) {
        let point = seeds[triangle[k]], {edges, vertices} = this.get(point), segment;
        let hasVertex = vertices.some(p => p.x === center.x && p.y === center.y);
        if (coCirculars.length === 0) { // clear triangulation
          if (hasVertex) {continue triangulation;};
          side: for (let pointer = 0; pointer < 3; pointer = pointer + 2) {// pointer for each of the two mediatrices per angle
            segment = segments[(k + pointer) % 3];
            if (segment in edges) {edges[segment]++;} else {edges[segment] = 1;};};
          vertices.push(center);}
        else { // ambiguous case as more cocyclic points (> 3 seeds)
          side: for (let pointer = 0; pointer < 3; pointer = pointer + 2) {// pointer for each of the two mediatrices for a point
            segment = segments[(k + pointer) % 3];
            let hasEdge = segment in edges;
            if (hasVertex && hasEdge) {continue side;};
            let interceptors = mediatrices[(k + pointer) % 3].intersect(center, radius); // circle/mediatrice two intersections
            if (interceptors.every(p => [seeds[triangle[(k + pointer + 2) % 3]], ...coCirculars].some(e => e.to(p) < point.to(p)))) {continue side;};
            if (hasEdge) {edges[segment]++;} else {edges[segment] = 1;};
            if (!hasVertex) {vertices.push(center);};};};};
    }; return this;}
  areas() {// Gauss formula (shoelace triangle form)
    const angle = (a, b, c) => { // angle ABC of vectors AB and BC (B is common point)
      let x = (a.x - b.x) * (c.x - b.x) + (a.y - b.y) * (c.y - b.y); // dot product
      let y = (a.x - b.x) * (c.y - b.y) + (c.x - b.x) * (a.y - b.y); // cross product
      return Math.atan2(y, x);}; // atan2 corrects error of angle versus atan(y/x) whenever x<0 (diametrical opposite)
    let areas = [];
    for (let [seed, {vertices, edges}] of this[Symbol.iterator]()) {
      let isOpenCell = Object.values(edges).some(ends => ends !== 2);
      let angles = vertices.length;
      if (isOpenCell || angles < 3) {areas.push(-1); continue;};
      let north = new Point(seed.x, seed.y + 1);
      vertices.sort((p,q) => (angle(q, seed, north) - angle(p, seed, north))); // sort counterclock wise from south
      let area = 0, i = 0;
      do {area += (vertices[i].x * vertices[(i+1) % angles].y) - (vertices[i].y * vertices[(i+1) % angles].x);} // product of determinants
      while (++i < angles);
      areas.push(area / 2);
    };
    return areas;}
};


// First Test
let squareCell = [new Point(0,0), new Point(2,0), new Point(-2,0), new Point(0,2), new Point(0,-2)];
// Second Test
let triangles = [new Point(2,1), new Point(2,-1), new Point(4.4,2.2),
  new Point(4.4,-2.2), new Point(-.4, 2.2), new Point(-.4, -2.2)];
// Last Test
let colinears = [new Point (1, 2), new Point (1.5, 3), new Point (2, 4), new Point (5, 10)];
// Square Grid Test
let squareGrid = (rows, cols) => {let output = [];
  for (let i = 0; i < rows; i++) {for (let j = 0; j < cols; j++) {output.push(new Point(j, i));};};
  return output;};
// Hexagonal Grid Test
let hexagonalGrid = (rows, cols) => {let output = [];
  for (let i = 0; i < rows; i++) {for (let j = 0; j < cols; j++) {output.push(new Point(i*1.7320508075688772, j*2+i));};};
  return output;};


let graph = new Voronoï(triangles); console.log(graph.cells().areas());
for (let [seed, {vertices, edges}] of graph[Symbol.iterator]()) {console.log(edges,/* ...vertices*/);};


/* To dos:
X - setting borders
X - computing mediatrices equations
X - computing triangulations
X - graph clean up for each seed point
O - convexHull => sorting points around centroïds
X - clearing ambiguous cases for triangulation of more than three cocyclic points
X - polyArea => shoelace computation
X - return areas or infinite mark (-1) */