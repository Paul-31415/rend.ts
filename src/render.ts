import * as PIXI from "pixi.js";
import 'pixi-display';
//import * as GPU from "gpu.js";

/*

optimizations:
[ ] - gpujs
[-] - caching:
  [√] - rays
  [√] - intersections
  [ ] - other
[ ] - colorer ray quantity
[ ] - other...

*/



import { RGBColor, Color } from "./color";
import { ucs2 } from "punycode";


const EPSILON = 0.00000001;

class V3D {
    constructor(public x: number,
        public y: number,
        public z: number) { }
    public times(k: number): V3D {
        return new V3D(this.x * k, this.y * k, this.z * k);
    }
    public minus(v: V3D): V3D {
        return new V3D(this.x - v.x, this.y - v.y, this.z - v.z);
    }
    public plus(v: V3D): V3D {
        return new V3D(this.x + v.x, this.y + v.y, this.z + v.z);
    }
    public dot(v: V3D): number {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }
    public cross(v: V3D): V3D {
        return new V3D(this.y * v.z - this.z * v.y, this.z * v.x - this.x * v.z, this.x * v.y - this.y * v.x);
    }
    public mag2(): number {
        return this.x * this.x + this.y * this.y + this.z * this.z;
    }
    public mag(): number {
        return Math.sqrt(this.mag2());
    }
    public norm(): V3D {
        var mag = this.mag();
        var div = (mag === 0) ? Infinity : 1. / mag;
        return this.times(div);
    }

    public max(v: V3D): V3D {
        return new V3D(Math.max(v.x, this.x), Math.max(v.y, this.y), Math.max(v.z, this.z));
    }
    public min(v: V3D): V3D {
        return new V3D(Math.min(v.x, this.x), Math.min(v.y, this.y), Math.min(v.z, this.z));
    }

    public static randSphere(): V3D {
        //adapted from : https://stackoverflow.com/questions/5531827/random-point-on-a-given-sphere
        var u: number = Math.random();
        var v: number = Math.random();
        var theta: number = 2 * Math.PI * u;
        var phi: number = Math.acos(2 * v - 1);
        var x: number = (Math.sin(phi) * Math.cos(theta));
        var y: number = (Math.sin(phi) * Math.sin(theta));
        var z: number = (Math.cos(phi));
        return new V3D(x, y, z);
    }

}

class Ray {
    constructor(public pos: V3D,
        public dir: V3D) { }
    public static intersectsBoundingBox(r: Ray, vmin: V3D, vmax: V3D): boolean {
        //adapted from https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms

        var t1: number = (vmin.x - r.pos.x) / r.dir.x;
        var t2: number = (vmax.x - r.pos.x) / r.dir.x;
        var t3: number = (vmin.y - r.pos.y) / r.dir.y;
        var t4: number = (vmax.y - r.pos.y) / r.dir.y;
        var t5: number = (vmin.z - r.pos.z) / r.dir.z;
        var t6: number = (vmax.z - r.pos.z) / r.dir.z;

        var tmin: number = Math.max(Math.max(Math.min(t1, t2), Math.min(t3, t4)), Math.min(t5, t6));
        var tmax: number = Math.min(Math.min(Math.max(t1, t2), Math.max(t3, t4)), Math.max(t5, t6));

        // if tmax < 0, ray (line) is intersecting AABB, but the whole AABB is behind us
        if (tmax < 0) {
            //t = tmax;
            return false;
        }

        // if tmin > tmax, ray doesn't intersect AABB
        if (tmin > tmax) {
            //t = tmax;
            return false;
        }

        //t = tmin;
        return true;

    }
}





/* abstractions:

Ray: class
  a ray in 3 space, no more

Traceable: interface
  traceables are anything that can be intersected by a ray

Intersection: class
  intersections contain information from tracing a ray to a traceable:
  they contain:
    traceable in question,
    whether it intersects,
	how far along the ray it intersects,
	color information:
	  a colorer object

colorer: interface
  controls how to calculate the color value from an intersection
  has method:
    getColor(scene,	the intersection in question, steps)

*/



interface Scene extends ColoredTraceable {

}
class ArrayScene implements Scene {
    constructor(public contents: Traceable[], public colorer = new NoColorer, public hasSky = false, public skyDistance = Infinity) { }
    intersection(r: Ray): Intersection {
        var I: Intersection = { thing: this, colorer: this.colorer, dist: this.skyDistance, intersects: this.hasSky, ray: r, work: null };
        for (var i = 0; i < this.contents.length; i++) {
            var I2 = this.contents[i].intersection(r);
            if (I2.intersects && I2.dist < I.dist) {
                I = I2;
            }
        }
        return I;
    }
    uvCoords(i: Intersection): number[] {
        //spherical coords at infinity
        var n: V3D = i.ray.dir.norm();
        var u: number = Math.atan2(n.x, n.z) / (2 * Math.PI) + 0.5;
        var v: number = n.y / 2 + .5;
        return [u, v];
    }
    normal(i: Intersection): V3D {
        //negative direction
        return i.ray.dir.times(-1);
    }


}




interface Intersection {
    thing: Traceable;
    intersects: boolean;
    //inside: boolean; implied with sign of r.dir dot normal 
    ray: Ray;
    dist: number;//minimum positive distance along the ray that intersects, or, if no intersection, closest point or something like that (unused yet)
    work: any; // saved work for use in u-v coords or normal
    colorer?: Colorer;
}
function noIntersection(thing: Traceable, r: Ray): Intersection {
    return { thing: thing, intersects: false, ray: r, colorer: new NoColorer(), dist: Infinity, work: null };
}



interface Traceable {
    intersection(r: Ray): Intersection;
}

interface ColoredTraceable extends Traceable {
    colorer: Colorer;
    uvCoords(i: Intersection): number[];
    normal(i: Intersection): V3D;

    //maximalPoint(): V3D
    //minimalPoint(): V3D;
}
class ColoredTraceableDefaults {
    constructor(public color: Color) { }
    intersection(r: Ray): Intersection {
        return noIntersection(this, r);
    }
    uvCoords(_i: Intersection): number[] {
        return [0, 0];
    }
    normal(_i: Intersection): V3D {
        return new V3D(1, 0, 0);
    }
    maximalPoint(): V3D {
        return new V3D(Infinity, Infinity, Infinity);
    }
    minimalPoint(): V3D {
        return new V3D(-Infinity, -Infinity, -Infinity);
    }
}


class Triangle implements ColoredTraceable {
    constructor(public colorer: Colorer,
        public pos: V3D,
        public ix: V3D,
        public iy: V3D) { }
    public intersection(r: Ray): Intersection {
        var h: V3D = r.dir.cross(this.iy);
        var a: number = this.ix.dot(h);
        if (a === 0) {
            return noIntersection(this, r);
        }
        var f: number = 1. / a;
        var s: V3D = r.pos.minus(this.pos)
        var u: number = f * (s.dot(h));
        if (u < 0 || u > 1) {
            return noIntersection(this, r);
        }
        var q: V3D = s.cross(this.ix);
        var v: number = f * r.dir.dot(q);
        if (v < 0 || u + v > 1) {
            return noIntersection(this, r);
        }
        var t: number = f * this.iy.dot(q);
        return { thing: this, intersects: t > 0, ray: r, dist: t, work: [u, v], colorer: this.colorer };
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            return i.work;
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    normal(_i: Intersection): V3D {
        return this.ix.cross(this.iy);
    }
    maximalPoint(): V3D {
        return this.pos.max(this.pos.plus(this.ix)).max(this.pos.plus(this.iy));
    }
    minimalPoint(): V3D {
        return this.pos.min(this.pos.plus(this.ix)).min(this.pos.plus(this.iy));
    }



}
class Parallelogram implements ColoredTraceable {
    constructor(public colorer: Colorer,
        public pos: V3D,
        public ix: V3D,
        public iy: V3D) { }
    public intersection(r: Ray): Intersection {
        var h: V3D = r.dir.cross(this.iy);
        var a: number = this.ix.dot(h);
        if (a === 0) {
            return noIntersection(this, r);
        }
        var f: number = 1. / a;
        var s: V3D = r.pos.minus(this.pos)
        var u: number = f * (s.dot(h));
        if (u < 0 || u > 1) {
            return noIntersection(this, r);
        }
        var q: V3D = s.cross(this.ix);
        var v: number = f * r.dir.dot(q);
        if (v < 0 || v > 1) {
            return noIntersection(this, r);
        }
        var t: number = f * this.iy.dot(q);
        return { thing: this, intersects: t > 0, ray: r, dist: t, work: [u, v], colorer: this.colorer };
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            return i.work;
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    normal(_i: Intersection): V3D {
        return this.ix.cross(this.iy);
    }
    maximalPoint(): V3D {
        return this.pos.max(this.pos.plus(this.ix)).max(this.pos.plus(this.iy)).max(this.pos.plus(this.iy).plus(this.ix));
    }
    minimalPoint(): V3D {
        return this.pos.min(this.pos.plus(this.ix)).min(this.pos.plus(this.iy)).min(this.pos.plus(this.iy).plus(this.ix));
    }

}

class Plane implements ColoredTraceable {
    constructor(public colorer: Colorer,
        public pos: V3D,
        public ix: V3D,
        public iy: V3D) { }
    public intersection(r: Ray): Intersection {
        var h: V3D = r.dir.cross(this.iy);
        var a: number = this.ix.dot(h);
        if (a === 0) {
            return noIntersection(this, r);
        }
        var f: number = 1. / a;
        var s: V3D = r.pos.minus(this.pos)
        var u: number = f * (s.dot(h));
        var q: V3D = s.cross(this.ix);
        var v: number = f * r.dir.dot(q);
        var t: number = f * this.iy.dot(q);
        return { thing: this, intersects: t > 0, ray: r, dist: t, work: [u, v], colorer: this.colorer };
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            return i.work;
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    normal(_i: Intersection): V3D {
        return this.ix.cross(this.iy);
    }
    maximalPoint(): V3D {
        return new V3D(Infinity, Infinity, Infinity);
    }
    minimalPoint(): V3D {
        return new V3D(-Infinity, -Infinity, -Infinity);
    }
}


class Sphere implements ColoredTraceable {
    constructor(public colorer: Colorer,
        public pos: V3D, public r: number) { }
    public intersection(r: Ray): Intersection {
        var L: V3D = r.pos.minus(this.pos);
        var a: number = r.dir.mag2();
        var b: number = 2 * r.dir.dot(L);
        var c: number = L.mag2() - this.r * this.r;

        if (a === 0) {
            return noIntersection(this, r);
        }
        var x = -b / (2 * a);
        var dx2 = (b * b - 4 * a * c);
        if (dx2 < 0) {
            return noIntersection(this, r);
        }
        var dx = Math.sqrt(dx2) / (2 * a);

        if (x - dx > 0) {
            return { thing: this, intersects: true, ray: r, dist: x - dx, work: null, colorer: this.colorer };
        }
        if (x + dx > 0) {
            return { thing: this, intersects: true, ray: r, dist: x + dx, work: null, colorer: this.colorer };
        }
        return { thing: this, intersects: false, ray: r, dist: x, work: null };
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            var n: V3D = i.ray.pos.plus(i.ray.dir.times(i.dist)).minus(this.pos).times(1 / this.r);
            var u: number = Math.atan2(n.x, n.z) / (2 * Math.PI) + 0.5;
            var v: number = n.y / 2 + .5;
            return [u, v];
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    normal(i: Intersection): V3D {
        if (i.thing === this) {
            var n: V3D = i.ray.pos.plus(i.ray.dir.times(i.dist)).minus(this.pos).times(1 / this.r);
            return n;
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    maximalPoint(): V3D {
        return this.pos.plus(new V3D(this.r, this.r, this.r));
    }
    minimalPoint(): V3D {
        return this.pos.plus(new V3D(-this.r, -this.r, -this.r));
    }
}

interface Transform {
    forward(r: Ray): Ray;
    forward(r: V3D): V3D;
    inverse(r: Ray): Ray;
    inverse(r: V3D): V3D;
}
class TransformedObj implements ColoredTraceable {
    constructor(public colorer: Colorer, public thing: ColoredTraceable,
        public transform: Transform) { }

    intersection(r: Ray): Intersection {
        var newRay: Ray = this.transform.inverse(r);
        var res: Intersection = this.thing.intersection(newRay);
        return { thing: this, intersects: res.intersects, ray: r, dist: res.dist, work: res };
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            return this.thing.uvCoords(i.work);
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    normal(i: Intersection): V3D {
        if (i.thing === this) {
            return this.transform.inverse(this.thing.normal(i.work));
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }
    }
    /*maximalPoint(): V3D {
        var M: V3D = this.thing.maximalPoint();
        var m: V3D = this.thing.minimalPoint();
        var v: V3D = this.transform.forward(M);
        for (var i = 1; i < 8; i++) {
            v = v.max(this.transform.forward(new V3D((i & 1) ? m.x : M.x, (i & 2) ? m.y : M.y, (i & 4) ? m.z : M.z)));
        }
        return v;
    }
    minimalPoint(): V3D {
        var M: V3D = this.thing.maximalPoint();
        var m: V3D = this.thing.minimalPoint();
        var v: V3D = this.transform.forward(M);
        for (var i = 1; i < 8; i++) {
            v = v.min(this.transform.forward(new V3D((i & 1) ? m.x : M.x, (i & 2) ? m.y : M.y, (i & 4) ? m.z : M.z)));
        }
        return v;
    }*/
}


/*class Group implements ColoredTraceable {
    constructor(public color: Color,
        public contents: ColoredTraceable[]) { }
    maximalPoint(): V3D {
        var v: V3D = new V3D(-Infinity, -Infinity, -Infinity);
        for (var i = 0; i < this.contents.length; i++) {
            v = v.max(this.contents[i].maximalPoint());
        }
        return v;
    }
    minimalPoint(): V3D {
        var v: V3D = new V3D(Infinity, Infinity, Infinity);
        for (var i = 0; i < this.contents.length; i++) {
            v = v.min(this.contents[i].minimalPoint());
        }
        return v;
    }
    intersection(r: Ray): Intersection {
        if (Ray.intersectsBoundingBox(this.minimalPoint(), this.maximalPoint())) {
            var I: Intersection = { thing: this, intersects: false, dist: Infinity, ray: r, work: null };
            var maxIndex = 0;
            for (var i = 0; i < this.contents.length; i++) {
                var I2: Intersection = this.contents[i].intersection(r);
                if (I2.intersects && I2.dist < I.dist) {
                    I = I2;
                    maxIndex = i;
                }
            }
            if (!I.intersects) {
                return { thing: this, intersects: false, dist: Infinity, ray: r, work: null };
            } else {
                return { thing: this, intersects: true, dist: I.dist, ray: r, work: [I, maxIndex] };
            }
        }
        else {
            return { thing: this, intersects: false, dist: 0, ray: r, work: null }
        }
    }
    uvCoords(i: Intersection): number[] {
        if (i.thing === this) {
            if (!i.intersects) {
                return [];
            } else {
                return [i.work[1]].concat(this.contents[i.work[1]].uvCoords(i.work[0]));
            }
        } else {
            throw new RangeError("Thing is not same as intersection thing");
        }

    }
}*/
interface Colorer {
    getColor(s: Scene, i: Intersection, steps: number): Color;
}
class NoColorer implements Colorer {
    constructor() { }
    getColor(_s: Scene, _i: Intersection, _steps: number): Color {
        return new RGBColor(0, 0, 0);
    }
}
class FlatColorer implements Colorer {
    constructor(public color: Color) { }
    getColor(_s: Scene, _i: Intersection, _steps: number): Color {
        return this.color;
    }
}
class LambertianColorer implements Colorer {
    constructor(public color: Color) { }
    getColor(s: Scene, i: Intersection, steps: number): Color {
        if (steps <= 0) {
            return this.color;
        } else {
            var v: V3D = V3D.randSphere();
            if (v.dot((i.thing as ColoredTraceable).normal(i)) < 0) {
                v = v.times(-1);
            }
            var r: Ray = new Ray(i.ray.pos.plus(i.ray.dir.times(i.dist)), v);
            r.pos = r.pos.plus(r.dir.times(EPSILON));
            var i2: Intersection = s.intersection(r);
            if (i2.colorer === undefined) {
                return this.color;
            }
            return this.color.apply(i2.colorer.getColor(s, i2, steps - 1));
        }
    }
}
class ReflectiveColorer implements Colorer {
    constructor(public color: Color) { }
    getColor(s: Scene, i: Intersection, steps: number): Color {
        if (steps <= 0) {
            return this.color;
        } else {
            var v: V3D = i.ray.dir;
            var norm: V3D = (i.thing as ColoredTraceable).normal(i).norm();
            v = v.plus(norm.times(-2 * norm.dot(v)));

            var r: Ray = new Ray(i.ray.pos.plus(i.ray.dir.times(i.dist)), v);
            r.pos = r.pos.plus(r.dir.times(EPSILON));
            var i2: Intersection = s.intersection(r);
            if (i2.colorer === undefined) {
                return this.color;
            }
            return this.color.apply(i2.colorer.getColor(s, i2, steps - 1));
        }
    }
}
class uvCheckerFlatColorer implements Colorer {
    constructor(public color1: Color,
        public color2: Color,
        public su: number,
        public sv: number) { }
    getColor(_s: Scene, i: Intersection, _steps: number): Color {
        var c: number[] = (i.thing as ColoredTraceable).uvCoords(i);
        if (((Math.floor(c[0] * this.su) + Math.floor(c[1] * this.sv)) & 1) == 0) {
            return this.color1;
        } else { return this.color2; }
    }
}
function mod(x: number, m: number): number {
    return ((x % m) + m) % m;
}
/*class uvFlatImageMap implements Colorer {
    constructor(public image: PIXI.Texture,
        public scaleU = 1,
        public scaleV = 1,
        public offsetU = 0,
        public offsetV = 0,
    ) { }
    getColor(_s: Scene, i: Intersection, _steps: number): Color {
        var c: number[] = (i.thing as ColoredTraceable).uvCoords(i);
        mod(c[0] * this.scaleU + this.offsetU, this.image.width);
        mod(c[1] * this.scaleV + this.offsetV, this.image.height);



    }

}*/





class uvCheckerColorer implements Colorer {
    constructor(public colorer1: Colorer,
        public colorer2: Colorer,
        public su: number,
        public sv: number) { }
    getColor(_s: Scene, i: Intersection, _steps: number): Color {
        var c: number[] = (i.thing as ColoredTraceable).uvCoords(i);
        if (((Math.floor(c[0] * this.su) + Math.floor(c[1] * this.sv)) & 1) == 0) {
            return this.colorer1.getColor(_s, i, _steps);
        } else { return this.colorer2.getColor(_s, i, _steps); }
    }
}


interface Texture {
    getColorer(): Colorer
}
interface UVTexture extends Texture {
    getUVColorer(u: number, v: number): Colorer
}

class MandelbrotTexture implements UVTexture {
    static iters = 40;
    getColorer(): Colorer {
        return new FlatColorer(new RGBColor(0, 0, 0));
    }

    getUVColorer(u: number, v: number): Colorer {
        var zr: number = u;
        var zi: number = v;
        var zrsq: number = u * u;
        var zisq: number = v * v;
        var i: number;
        for (i = 0; zrsq + zisq < 4 && i < MandelbrotTexture.iters; i++) {
            zi = zr * zi * 2 + v;
            zr = zrsq - zisq + u;
            zrsq = zr * zr;
            zisq = zi * zi;
        }
        if (i < MandelbrotTexture.iters) {
            return new LambertianColorer(new RGBColor((Math.sin(i) + 1) / 3, (Math.sin(i * (2 ** .5)) + 1) / 3, (Math.sin(i * (3 ** .5)) + 1) / 3));
        } else {
            return new ReflectiveColorer(new RGBColor(.9, .9, .9));
        }
    }
}


class uvColorer implements Colorer {
    constructor(public texture: UVTexture,
        public scaleU = 1,
        public scaleV = 1,
        public offsetU = 0,
        public offsetV = 0,
    ) { }
    getColor(_s: Scene, i: Intersection, _steps: number): Color {
        var c: number[] = (i.thing as ColoredTraceable).uvCoords(i);
        return this.texture.getUVColorer(c[0] * this.scaleU + this.offsetU, c[1] * this.scaleV + this.offsetV).getColor(_s, i, _steps);
    }

}



class ColorerAverage implements Colorer {
    constructor(public colorer1: Colorer, public w1: number,
        public colorer2: Colorer) { }
    getColor(_s: Scene, i: Intersection, _steps: number): Color {
        if (Math.random() < this.w1) {
            return this.colorer1.getColor(_s, i, _steps);
        } else {
            return this.colorer2.getColor(_s, i, _steps);
        }
    }
}


function traceOrtho(scene: Scene, start: Ray, ix: V3D, iy: V3D, dx: number, dy: number, steps: number = 10, samples: number = 1000): PIXI.Graphics {
    var result: PIXI.Graphics = new PIXI.Graphics();
    var dix: V3D = ix.times(dx);
    var diy: V3D = iy.times(dy);

    //result.addChild(new PIXI.Text("test", { fill: 0xffffff }));

    var r: Ray = start;
    var y: number = 0;
    for (var fy = 0; fy < 1; fy += dy) {
        var R: Ray = new Ray(r.pos, r.dir);
        result.moveTo(0, y);
        var x: number = 0;
        for (var fx = 0; fx < 1; fx += dx) {
            var i: Intersection = scene.intersection(R);


            var color: number = ((x + y) & 1) ? 0x808080 : 0xff00ff;


            if (i.colorer != undefined) {
                var c: Color = i.colorer.getColor(scene, i, steps);
                for (var s = 1; s < samples; s++) {
                    c = c.plus(i.colorer.getColor(scene, i, steps));
                }
                color = c.scale(1. / samples).getRGB();

            }

            result.lineStyle(1, color, 1, 1);
            if (i.intersects) {
                result.lineTo(x, y);
            } else {
                result.moveTo(x, y);
            }

            R.pos = R.pos.plus(dix);
            x++;
        }
        r.pos = r.pos.plus(diy);
        y++;
    }
    return result;
}

function tracePersp(scene: Scene, start: Ray, ix: V3D, iy: V3D, dx: number, dy: number, steps: number = 10, samples: number = 1000): ColorImage {
    var result: ColorImage = new ColorImage(Math.ceil(1 / dx), Math.ceil(1 / dy));
    var dix: V3D = ix.times(dx);
    var diy: V3D = iy.times(dy);

    //result.addChild(new PIXI.Text("test", { fill: 0xffffff }));

    var r: Ray = start;
    var y: number = 0;
    for (var fy = 0; fy < 1; fy += dy) {
        var R: Ray = new Ray(r.pos, r.dir);
        var x: number = 0;
        for (var fx = 0; fx < 1; fx += dx) {
            var i: Intersection = scene.intersection(R);


            var color: Color = ((x + y) & 1) ? new RGBColor(.5, .5, .5) : new RGBColor(1, 0, 1);


            if (i.colorer != undefined) {
                var c: Color = i.colorer.getColor(scene, i, steps);
                for (var s = 1; s < samples; s++) {
                    c = c.plus(i.colorer.getColor(scene, i, steps));
                }
                color = c;

            }

            result.colors[x][y] = new SampledColor(color.scale(1 / samples), samples, R);

            R.dir = R.dir.plus(dix);
            x++;
        }
        r.dir = r.dir.plus(diy);
        y++;
    }
    return result;
}

function traceRefine(scene: Scene, image: ColorImage, steps: number = 10, sampleFactor: number = 10) {
    for (var y = 0; y < image.height; y++) {
        for (var x = 0; x < image.width; x++) {
            var i: Intersection = scene.intersection(image.colors[x][y].r);


            var color: Color = ((x + y) & 1) ? new RGBColor(.5, .5, .5) : new RGBColor(1, 0, 1);


            if (i.colorer != undefined) {
                var samples: number = 1;
                var c: Color = i.colorer.getColor(scene, i, steps);
                for (var s = 1; s < sampleFactor; s += 1 + 0 * image.colors[x][y].c.brightness()) {
                    c = c.plus(i.colorer.getColor(scene, i, steps));
                    samples++;
                }
                color = c;
                image.colors[x][y].addColor(color, samples);
            }




        }
    }
}
/*
var gpu = new GPU.GPU();
function traceRefineKernal(s: Scene, I: ColorImage, steps = 10) {
    var i: Intersection = s.intersection(I.colors[this.thread.x][this.thread.y].r);
    if (i.colorer != undefined) {
        I.colors[this.thread.x][this.thread.y].addColor(i.colorer.getColor(s, i, steps), 1);
    }
}
const traceRefineGPU = gpu.createKernel(traceRefineKernal);
*/


var colr1: Colorer = new FlatColorer(new RGBColor(30, 25, 20));
var colr2: Colorer = new FlatColorer(new RGBColor(0, 0.5, 200));
var colr3: Colorer = new ColorerAverage(new ReflectiveColorer(new RGBColor(0, 0.1, 0)), 0.9, new LambertianColorer(new RGBColor(0, 0, 1)));


var light: Sphere = new Sphere(colr1, new V3D(20, -20, 20), 6);
var ground: Plane = new Plane(/*new uvCheckerColorer(new ColorerAverage(new ReflectiveColorer(new RGBColor(0.1, 0.1, 0.1)), 0.9, new LambertianColorer(new RGBColor(0.1, 0.1, .1))), new ColorerAverage(new ReflectiveColorer(new RGBColor(0.8, 0.8, 0.8)), 0.9, new LambertianColorer(new RGBColor(0.8, 0.8, .8))), 1, 1)*/ new uvColorer(new MandelbrotTexture(), 0.1, 0.1, -.4, -.6), new V3D(0, 2, 0), new V3D(1, 0, 1), new V3D(-1, 0, 1));
var sphere1: Sphere = new Sphere(colr3, new V3D(0, 0, 0), 0.5);
var sphere2: Sphere = new Sphere(colr3, new V3D(0, 1.2, 4), 2);
var triangle1: Triangle = new Triangle(colr2, new V3D(0, 0, 4), new V3D(-1, -1, -1), new V3D(-1, 0, 0));
var triangle2: Triangle = new Triangle(new FlatColorer(new RGBColor(100, 0.5, 0)), new V3D(-20, -20, 100), new V3D(-110, -13, 120), new V3D(-133, 30, 30));

var exampleScene: Scene = new ArrayScene([sphere1, triangle1, triangle2, sphere2, ground, light], new uvCheckerFlatColorer(new RGBColor(0.125, 0.125, 0.125), new RGBColor(0, 0, 0), 10, 10), true);

for (var y = 2; y > -5; y -= 0.1) {
    (exampleScene as ArrayScene).contents = (exampleScene as ArrayScene).contents.concat([
        new Triangle(colr3, new V3D(2, y, 2), new V3D(1, 0.1, 0), new V3D(0, 0.1, 1))]
    );
}






class SampledColor {
    public r: Ray;
    constructor(public c: Color, public s: number, r: Ray) {
        this.r = new Ray(new V3D(r.pos.x, r.pos.y, r.pos.z), new V3D(r.dir.x, r.dir.y, r.dir.z));
    }
    public addColor(c: Color, s: number) {
        this.c = this.c.scale(this.s).plus(c).scale(1 / (this.s + s));
        this.s += s;
    }
    public add(c: SampledColor) {
        this.c = this.c.scale(this.s).plus(c.c.scale(c.s)).scale(1 / (this.s + c.s));
        this.s += c.s;
    }
}


class ColorImage {
    public colors: SampledColor[][];

    constructor(public width: number, public height: number) {
        this.colors = new Array(width)
        for (var i = 0; i < width; i++) {
            this.colors[i] = new Array(height);
        }

    }

    public add(o: ColorImage) {
        for (var y = 0; y < this.height; y++) {
            for (var x = 0; x < this.width; x++) {
                this.colors[x][y].add(o.colors[x][y]);
            }
        }

    }


    public get(scale: number): PIXI.Graphics {
        var result: PIXI.Graphics = new PIXI.Graphics();
        for (var y = 0; y < this.height; y++) {
            result.moveTo(0, y * scale);
            for (var x = 0; x < this.width; x++) {
                result.lineStyle(scale, this.colors[x][y].c.getRGB(),
                    this.colors[x][y].c.getAlpha(), 1);
                result.lineTo(x * scale, y * scale);
            }
        }
        return result;
    }
}




export {
    Ray, V3D, ColorImage,
    traceOrtho, tracePersp, traceRefine,

    Traceable,
    Intersection,
    ColoredTraceable,

    Triangle, Parallelogram, Plane, Sphere,

    Texture, UVTexture,

    Colorer,
    NoColorer, FlatColorer, LambertianColorer, ReflectiveColorer, uvColorer,

    exampleScene
}
