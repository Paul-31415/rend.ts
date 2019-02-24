import * as PIXI from "pixi.js";



interface Angle {
    apply(x: number, y: number, z: number): number[];
    plus(o: Angle): Angle;
    minus(o: Angle): Angle;
    inverse(): Angle;
    scale(n: number): Angle;
    plusV(v: number[], dt: number): Angle;
}
/*
class EulerAngle {
    constructor(public x = 0, public y = 0, public z = 0) { }
    scale(n = 1): Angle {
        return new EulerAngle(n * this.x, n * this.y, n * this.z);
    }
    inverse(): Angle {
        return this.scale(-1);
    }
    minus(o: Angle): Angle {
        return this.plus(o.inverse());
    }
    dot(x = 0, y = 0, z = 0): number {
        return this.x * x + this.y * y + this.z * z;
    }
    add(o: EulerAngle): EulerAngle {
        return new EulerAngle(this.x + o.x, this.y + o.y, this.z + o.z);
    }
    cross(x = 0, y = 0, z = 0): number[] {
        return [
            this.y * z - y * this.z,
            this.z * x - z * this.x,
            this.x * y - x * this.y
        ];
    }
    mag(): number {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }
    apply(x = 0, y = 0, z = 0): number[] {
        const m = this.mag();
        if (m == 0) {
            return [x, y, z];
        }
        const c = Math.cos(m);
        const s = Math.sin(m);
        const e = this.scale(1 / m) as EulerAngle;
        const crs = e.cross(x, y, z);
        const dot = e.dot(x, y, z);
        return [c * x + s * crs[0] + (1 - c) * dot * e.x, c * x + s * crs[1] + (1 - c) * dot * e.y, c * x + s * crs[2] + (1 - c) * dot * e.z];
    }
    plus(a: Angle): Angle {
        const o = a as EulerAngle;
        const m = this.mag();
        const mc = Math.cos(m / 2);
        const ms = Math.sin(m / 2);
        const om = o.mag();
        const omc = Math.cos(om / 2);
        const oms = Math.sin(om / 2);
        const n = this.scale(1 / m) as EulerAngle;
        const on = o.scale(1 / om) as EulerAngle;


        const rots = 2 * Math.acos(mc * omc - ms * oms * n.dot(on.x, on.y, on.z));
        const c = on.cross(n.x, n.y, n.z);
        const f = rots / Math.sin(rots / 2);
        return new EulerAngle((on.x * oms * mc + n.x * ms * omc + ms * oms * c[0]) * f, (on.y * oms * mc + n.y * ms * omc + ms * oms * c[1]) * f, (on.z * oms * mc + n.z * ms * omc + ms * oms * c[2]) * f);
    }
}// */

//*
class QuaternionAngle {
    constructor(public r = 1, public i = 0, public j = 0, public k = 0) { }
    times(o: QuaternionAngle): QuaternionAngle {
        return new QuaternionAngle(
            this.r * o.r - this.i * o.i - this.j * o.j - this.k * o.k,
            this.r * o.i + this.i * o.r + this.j * o.k - this.k * o.j,
            this.r * o.j + this.j * o.r + this.k * o.i - this.i * o.k,
            this.r * o.k + this.k * o.r + this.i * o.j - this.j * o.i);
    }
    dot(o: QuaternionAngle): number {
        return this.r * o.r + this.i * o.i + this.j * o.j + this.k * o.k;
    }
    add(o: QuaternionAngle): QuaternionAngle {
        return new QuaternionAngle(this.r + o.r, this.i + o.i, this.j + o.j, this.k + o.k);
    }
    static fromAxisAngle(x = 0, y = 0, z = 0): QuaternionAngle {
        const mag = Math.sqrt(x * x + y * y + z * z);
        if (mag < 0.00000000001) {
            return new QuaternionAngle(1, 0, 0, 0);
        }
        const f = Math.sin(mag / 2) / mag;
        return new QuaternionAngle(Math.cos(mag / 2), f * x, f * y, f * z);
    }
    plus(o: Angle): Angle {
        return (o as QuaternionAngle).times(this);
    }
    plusV(v: number[], dt: number): Angle {
        const mag = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * dt;
        if (mag < 0.000000001) {
            return this;
        }
        const f = Math.sin(mag / 2) / mag;
        return this.plus(new QuaternionAngle(Math.cos(mag / 2), v[0] * f, v[1] * f, v[2] * f));
        //return this.add(new QuaternionAngle(0, v[0] * dt / 2, v[1] * dt / 2, v[2] * dt / 2).times(this)).normalized();
    }

    inverse(): Angle {
        const mag2 = this.r * this.r + this.i * this.i + this.j * this.j + this.k * this.k;
        return new QuaternionAngle(this.r / mag2, -this.i / mag2, -this.j / mag2, -this.k / mag2);
    }
    conjugate(): QuaternionAngle {
        return new QuaternionAngle(this.r, -this.i, -this.j, -this.k);
    }
    apply(x = 0, y = 0, z = 0): number[] {
        const res = this.times(new QuaternionAngle(0, x, y, z)).times((this.inverse()) as QuaternionAngle);
        return [res.i, res.j, res.k];

    }
    minus(o: Angle): Angle {
        return this.plus((o as QuaternionAngle).conjugate());
    }
    scale(n = 1): Angle { //to be used only for angular vel
        return new QuaternionAngle(this.r * n, this.i * n, this.j * n, this.k * n);

    }
    normalized(): QuaternionAngle {
        const mag = Math.sqrt(this.r * this.r + this.i * this.i + this.j * this.j + this.k * this.k);
        return new QuaternionAngle(this.r / mag, this.i / mag, this.j / mag, this.k / mag);
    }
    axis(): number[] {
        const mag = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
        if (mag == 0) {
            return [1, 0, 0];
        } else {
            return [this.i / mag, this.j / mag, this.k / mag];
        }
    }
    angle(): number {
        const mag = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
        return 2 * Math.atan2(mag, this.r);
    }
    axisAngle(): number[] {
        const mag = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
        const theta = 2 * Math.atan2(mag, this.r);
        const s = Math.sin(theta / 2);
        if (s < 0.00000000001) {
            return [0, 0, 0];
        }
        /*if (Math.abs(theta) > Math.PI) {
            const f = (theta + (theta > 0 ? -2 : 2) * Math.PI) / theta;
            return [f * this.i / s, f * this.j / s, f * this.k / s];
        }*/
        return [theta * this.i / s, theta * this.j / s, theta * this.k / s];
    }

}
//debugger;
// */
/*
class MatrixAngle implements Angle {
    constructor(public x = [1, 0, 0], public y = [0, 0, 1], public z = [0, 0, 1]) { }
    apply(x = 0, y = 0, z = 0): number[] {
        return [x * this.x[0] + y * this.y[0] + z * this.z[0],
        x * this.x[1] + y * this.y[1] + z * this.z[1],
        x * this.x[2] + y * this.y[2] + z * this.z[2]];
    }
    inverse(): Angle {
        return new MatrixAngle([this.x[0], this.y[0], this.z[0]], [this.x[1], this.y[1], this.z[1]], [this.x[2], this.y[2], this.z[2]]);
    }
    plus(o: Angle): Angle {
        o = (o as MatrixAngle);
        return new MatrixAngle(o.apply(this.x), o.apply(this.y), o.apply(this.z));
    }
    minus(o: Angle): Angle {
        o = (o as MatrixAngle);
        return this.plus(o.inverse());
    }
}*/
/*    minus(o: Angle): Angle {
        return new Angle(this.r - o.r, this.i - o.i, this.j - o.j, this.k - o.k);
    }
    constructor(public r = 0, public i = 1, public j = 0, public k = 0) { }
    public static fromEuler(ex = 0, ey = 0, ez = 0, theta = 0): Angle {
        const m = Math.sqrt(ex * ex + ey * ey + ez * ez);
        const s = Math.sin(theta / 2);
        return new Angle(Math.cos(theta / 2), ex * s / m, ey * s / m, ez * s / m);
    }
    public static fromVector(ex = 0, ey = 0, ez = 0): Angle {
        const m = Math.sqrt(ex * ex + ey * ey + ez * ez);
        const s = Math.sin(m / 2);
        return new Angle(Math.cos(m / 2), ex * s / m, ey * s / m, ez * s / m);
    }
    public static AVfromVector(ex = 0, ey = 0, ez = 0): Angle {
        return new Angle(0, ex, ey, ez);
    }
 
    plus(o: Angle): Angle {
        return new Angle(this.r + o.r, this.i + o.i, this.j + o.j, this.k + o.k);
    }
    times(o: Angle): Angle {
        return new Angle(
            this.r * o.r - this.i * o.i - this.j * o.j - this.k * o.k,
            this.r * o.i + this.i * o.r + this.j * o.k - this.k * o.j,
            this.r * o.j + this.j * o.r + this.k * o.i - this.i * o.k,
            this.r * o.k + this.k * o.r + this.i * o.j - this.j * o.i);
    }
    conj(): Angle {
        return new Angle(this.r, -this.i, -this.j, -this.k);
    }
    conjugate(o: Angle): Angle {
        return this.times(o).times(this.reciprocal());
    }
    mag2(): number {
        return this.r * this.r + this.i * this.i + this.j * this.j + this.k * this.k;
    }
    scale(n: number): Angle {
        return new Angle(this.r * n, this.i * n, this.j * n, this.k * n);
    }
    inverse(): Angle {
        return new Angle(-this.r, this.i, this.j, this.k);
    }
    reciprocal(): Angle {
        return this.scale(1 / this.mag2());
    }
    mag(): number {
        return Math.sqrt(this.mag2());
    }
    unit(): Angle {
        return this.scale(1 / this.mag());
    }
    rotate(x = 0, y = 0, z = 0): number[] {
        const r = this.conjugate(new Angle(0, x, y, z));
        return [r.i, r.j, r.k];
    }
    toEuler(): number[] {
        const m = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
        return [this.i / m, this.j / m, this.k / m, 2 * Math.atan2(m, this.r)];
    }
    compose(o: Angle): Angle {
        return o.times(this);
    }
    applyVelocity(o: Angle, dt = 0): Angle {
        return this.plus(o.scale(dt / 2).times(this));
    }
*/










const defaultAngle: Angle = new QuaternionAngle();


class PPoint3d extends PIXI.Point {
    constructor(x: number, y: number, public z: number, public p: Physics, public a = defaultAngle) {
        super(x, y);
    }
    update(dt: number): void {
        this.p.update(this, dt);
    }
}

interface Physics {
    update(p: PPoint3d, dt: number): void;
    addImpulse(i: Impulse): void;
    vx: number; vy: number; vz: number;
    va: number[];
}

class noPhysics implements Physics {
    public vx = 0; public vy = 0; public vz = 0;
    public va = [0, 0, 0];
    update(_p: PPoint3d, _dt: number): void { }
    addImpulse(_t: Impulse): void { };
}
class Hinge implements Physics {
    public vx = 0; public vy = 0; public vz = 0;
    public va = [0, 0, 0];
    constructor(public I = 0) { }
    update(p: PPoint3d, dt: number): void {
        p.a = p.a.plusV(this.va, dt);
    }

    addImpulse(i: Impulse): void {
        this.va[0] += i.iax / this.I;
        this.va[1] += i.iay / this.I;
        this.va[2] += i.iaz / this.I;
    }
}

class GravityDecorator implements Physics {
    constructor(public gx: number, public gy: number, public gz: number, public p: Physics) { }
    update(P: PPoint3d, dt: number): void {
        this.p.update(P, dt);
        this.p.vx += this.gx * dt;
        this.p.vy += this.gy * dt;
        this.p.vz += this.gz * dt;
    }
    addImpulse(i: Impulse): void { this.p.addImpulse(i); }
    get va(): number[] { return this.p.va; }
    get vx() { return this.p.vx; }
    get vy() { return this.p.vy; }
    get vz() { return this.p.vz; }
    set va(n: number[]) { this.p.va = n; }
    set vx(n: number) { this.p.vx = n; }
    set vy(n: number) { this.p.vy = n; }
    set vz(n: number) { this.p.vz = n; }
}


class PVel implements Physics {
    constructor(public vx = 0, public vy = 0, public vz = 0, public va = [0, 0, 0]) { }
    update(p: PPoint3d, dt: number): void {
        p.x += this.vx * dt;
        p.y += this.vy * dt;
        p.z += this.vz * dt;
        p.a = p.a.plusV(this.va, dt);
    }
    addImpulse(_t: Impulse): void { };
}

interface Force {
    doForce(p: PPoint3d, dt: number): void;
}


class PNeut extends PVel {
    constructor(vx = 0, vy = 0, vz = 0, public mass = 1, public f: Force[] = [], va = [0, 0, 0], public I = 1) { super(vx, vy, vz, va); }
    update(p: PPoint3d, dt: number): void {
        for (var i = 0; i < this.f.length; i++) {
            this.f[i].doForce(p, dt);
        }
        p.x += this.vx * dt;
        p.y += this.vy * dt;
        p.z += this.vz * dt;
        p.a = p.a.plusV(this.va, dt);
    }

    addImpulse(i: Impulse): void {
        this.vx += i.ix / this.mass;
        this.vy += i.iy / this.mass;
        this.vz += i.iz / this.mass;
        this.va[0] += i.iax / this.I;
        this.va[1] += i.iay / this.I;
        this.va[2] += i.iaz / this.I;
    }
}
function dot(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number): number { return x1 * x2 + y1 * y2 + z1 * z2; }

class PLinkage extends PNeut {
    constructor(public link: PPoint3d, public l = 1, public stiffness = 1, public dampening = 0, vx = 0, vy = 0, vz = 0, mass = 1) { super(vx, vy, vz, mass) }
    symmetricAddForce(fx: number, fy: number, fz: number) {
        this.link.p.addImpulse(new Impulse(-fx, -fy, -fz));
        this.addImpulse(new Impulse(fx, fy, fz));
    }
    update(p: PPoint3d, dt: number): void {
        const dx = p.x - this.link.x;
        const dy = p.y - this.link.y;
        const dz = p.z - this.link.z;
        const dvx = this.vx - (this.link.p as PVel).vx;
        const dvy = this.vy - (this.link.p as PVel).vy;
        const dvz = this.vz - (this.link.p as PVel).vz;
        const dvDotDir = dot(dx, dy, dz, dvx, dvy, dvz) / this.l;
        const l2 = this.l * this.l;
        const len2 = dx * dx + dy * dy + dz * dz;
        const dl2 = l2 - len2;
        const fmag1 = (this.stiffness * dl2) * dt;
        const ix = dx / this.l;
        const iy = dy / this.l;
        const iz = dz / this.l;
        const fmag2 = - (this.dampening * dvDotDir) * dt;
        this.symmetricAddForce(ix * fmag1, iy * fmag1, iz * fmag1);
        p.x += this.vx * dt;
        p.y += this.vy * dt;
        p.z += this.vz * dt;
        this.symmetricAddForce(ix * fmag2, iy * fmag2, iz * fmag2);

    }
}


class FLinkage implements Force {
    symmetricAddForce(p: Physics, fx: number, fy: number, fz: number) {
        this.link.p.addImpulse(new Impulse(-fx, -fy, -fz));
        p.addImpulse(new Impulse(fx, fy, fz));
    }
    doForce(p: PPoint3d, dt: number): void {
        const dx = p.x - this.link.x;
        const dy = p.y - this.link.y;
        const dz = p.z - this.link.z;
        const dvx = p.p.vx - (this.link.p as PVel).vx;
        const dvy = p.p.vy - (this.link.p as PVel).vy;
        const dvz = p.p.vz - (this.link.p as PVel).vz;
        const dvDotDir = dot(dx, dy, dz, dvx, dvy, dvz) / this.l;
        const l2 = this.l * this.l;
        const len2 = dx * dx + dy * dy + dz * dz;
        const dl2 = l2 - len2;
        const fmag1 = (this.stiffness * dl2) * dt;
        const ix = dx / this.l;
        const iy = dy / this.l;
        const iz = dz / this.l;
        const fmag2 = - (this.dampening * dvDotDir) * dt;
        this.symmetricAddForce(p.p, ix * fmag1, iy * fmag1, iz * fmag1);
        //p.x += this.vx * dt;
        //p.y += this.vy * dt;
        //p.z += this.vz * dt;
        this.symmetricAddForce(p.p, ix * fmag2, iy * fmag2, iz * fmag2);

    }
    constructor(public link: PPoint3d, public l = 1, public stiffness = 1, public dampening = 0) { }

}

function cross(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number): number[] {
    return [
        y1 * z2 - y2 * z1,
        z1 * x2 - z2 * x1,
        x1 * y2 - x2 * y1
    ];


}
/*
class Vector{
	constructor(public x=0,public y=0,public z=0){}
	cross(o:Vector):Vector{
		return new Vector(   
        this.y * o.z - o.y * this.z,
        this.z * o.x - o.z * this.x,
        this.x * o.y - o.x * this.y
		);
	}
	dot(o:Vector):number{
		return this.x*o.x+this.y*o.y+this.z*o.z;
	}
	plus(o:Vector):Vector{
		return new Vector(this.x+o.x,this.y+ o.y, this.z+o.z);
	}
	minus(o:Vector):Vector{
		return new Vector(this.x-o.x,this.y- o.y, this.z-o.z);
	}
	times(o:Vector):Vector{
		return new Vector(this.x*o.x,this.y* o.y, this.z*o.z);
	}
	norm
}*/

class Impulse {
    constructor(public ix = 0, public iy = 0, public iz = 0,
        public iax = 0, public iay = 0, public iaz = 0) { }
    negative(): Impulse {
        return new Impulse(-this.ix, -this.iy, -this.iz, -this.iax, -this.iay, -this.iaz);
    }
    shift(x = 0, y = 0, z = 0): Impulse {//calculates the equivilent impulse for if this happens at x,y,z from the origin
        // only changes torque
        const it = cross(x, y, z, this.ix, this.iy, this.iz); // torque
        return new Impulse(this.ix, this.iy, this.iz, it[0] + this.iax, it[1] + this.iay, it[2] + this.iaz);
    }
    scale(n = 1): Impulse {
        return new Impulse(this.ix * n, this.iy * n, this.iz * n, this.iax * n, this.iay * n, this.iaz * n);
    }
    plus(i: Impulse): Impulse {//add two impulses at the same place
        return new Impulse(this.ix + i.ix, this.iy + i.iy, this.iz + i.iz, this.iax + i.iax, this.iay + i.iay, this.iaz + i.iaz);

    }
}


function internalImpulse(p1: PPoint3d, p2: PPoint3d, i: Impulse) {
    //conserve linear and angular momentum of the system






    //i is an impulse from p1 to p2
    //first, add the intended impulse to p2
    p2.p.addImpulse(i);
    //then add the shifted impulse to p1;
    p1.p.addImpulse(i.negative().shift(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z));
}



/*function pointAt(a: QuaternionAngle, x = 1, y = 0, z = 0): Angle {
    const m = Math.sqrt(x * x + y * y + z * z);
    const v = cross(a.i, a.j, a.k, x / m, y / m, z / m);
    const theta = Math.arccos(dot(a.i, a.j, a.k, x / m, y / m, z / m));

}*/

//debugger;
class FAngleLinkage implements Force {
    doForce(p: PPoint3d, dt: number): void {

        const ra1 = this.link.a.plus(this.a1);
        const ra2 = p.a.plus(this.a2);

        //go into where it's pointing

        //t = rxf
        //f = txr/(r^2)
        const dx = -(p.x - this.link.x);
        const dy = -(p.y - this.link.y);
        const dz = -(p.z - this.link.z);
        const r2 = dx * dx + dy * dy + dz * dz;
        const r = Math.sqrt(r2);

        // just find which way to push it
        const rp = ra1.apply(1, 0, 0);
        const tp = cross(dx / r, dy / r, dz / r, rp[0], rp[1], rp[2]);

        const ta = (ra2.minus(ra1) as QuaternionAngle).axisAngle();

        const ts = [tp[0] + ta[0], tp[1] + ta[1], tp[2] + ta[2]];


        //damping: attempt to match our angular velocity to link's
        // w = v/r 
        // v = rw

        const dvx = ((p.p as PVel).vx - (this.link.p as PVel).vx);
        const dvy = ((p.p as PVel).vy - (this.link.p as PVel).vy);
        const dvz = ((p.p as PVel).vz - (this.link.p as PVel).vz);
        const va = cross(dvx, dvy, dvz, dx / r2, dy / r2, dz / r2);

        const td = [this.link.p.va[0] - va[0], this.link.p.va[1] - va[1], this.link.p.va[2] - va[2]]

        const K = this.dampening * dt;

        const k = -this.stiffness * dt;

        const t = [ts[0] * k + td[0] * K, ts[1] * k + td[1] * K, ts[2] * k + td[2] * K];

        const f = cross(t[0], t[1], t[2], dx / r2, dy / r2, dz / r2);
        //t = rxf
        const et = cross(dx, dy, dz, f[0], f[1], f[2]);



        const IRes = new Impulse(f[0], f[1], f[2], t[0] - et[0], t[1] - et[1], t[2] - et[2]);


        internalImpulse(p, this.link, IRes);

    }
    constructor(public link: PPoint3d, public a1: Angle = defaultAngle, public a2: Angle = defaultAngle, public stiffness = 1, public dampening = 0) { }
}



/*
class CombinePhysics implements Physics {
    public p: Physics[];
    public vx: number;
    public vy: number;
    public vz: number;
    public va: Angle;

    constructor(...p: Physics[]) { this.p = p; this.vx = 0; this.vy = 0; this.vz = 0; this.va = defaultAngle; }
    update(p: PPoint3d, dt: number): void {
        for (var i = 0; i < this.p.length; i++) {
            this.p[i].update(p, dt);
        }
    }
    addForce(fx: number, fy: number, fz: number): void {
        for (var i = 0; i < this.p.length; i++) {
            this.p[i].addForce(fx, fy, fz);
        }
    }
    addTorque(t: Angle): void {
        for (var i = 0; i < this.p.length; i++) {
            this.p[i].addTorque(t);
        }
    }

}
*/

const g = 0;
function Tube(circ: number, links: number, len = 1, s = 1, d = 1, sc = 1, dc = 1, sups = 0, supd = 1, ): PPoint3d[] {
    var points: PPoint3d[] = [];
    const r = len / (2 * Math.sin(Math.PI / circ));
    for (var i = 0; i < circ; i++) {
        points.push(new PPoint3d(0, r * Math.sin(i * Math.PI * 2 / circ), r * Math.cos(i * Math.PI * 2 / circ), new noPhysics()));
    }


    const apothem = (len / 2) * Math.tan(Math.PI * (circ - 2) / (2 * circ));

    const k = Math.sqrt(((3 ** .5) / 2) ** 2 - ((r - apothem) / len) ** 2);

    /*
         5
      2     
         3
      0    
         4
      1    
         5
      2    
         3
 
         3:0,2
         4
         
      */
    for (var l = 1; l < links; l++) {
        points.push(new PPoint3d(l * len * k, r * Math.sin((0 - .5 * l) * Math.PI * 2 / circ), r * Math.cos((0 - .5 * l) * Math.PI * 2 / circ), new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - circ], len, sc, dc), new FLinkage(points[points.length - 1], len, sc, dc)]))))
        for (var i = 1; i < circ; i++) {
            points.push(new PPoint3d(l * len * k, r * Math.sin((i - .5 * l) * Math.PI * 2 / circ), r * Math.cos((i - .5 * l) * Math.PI * 2 / circ), new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - circ], len, sc, dc), new FLinkage(points[points.length - circ - 1], len, sc, dc), new FLinkage(points[points.length - 1], len, s, d)]))))
        }
        ((points[points.length - circ].p as GravityDecorator).p as PNeut).f.push(new FLinkage(points[points.length - 1], len, s, d))
        if (sups != 0) {
            for (var i = 0; i < circ / 2; i++) {
                ((points[points.length - i - (circ / 2) - 1].p as GravityDecorator).p as PNeut).f.push(new FLinkage(points[points.length - i - 1], 2 * r, sups, supd))
            }
        }

    }
    return points;
}




function TriBridge(links: number, len = 1, s = 1, d = 1, sc = 1, dc = 1): PPoint3d[] {
    const k = 3 ** .5 / 2;
    var points = [new PPoint3d(0, 0, 0, new noPhysics()), new PPoint3d(len / 2, len * k, 0, new noPhysics())];
    for (var i = 1; i < links; i++) {
        points.push(new PPoint3d(i * len, 0, 0, new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - 1], len, sc, dc), new FLinkage(points[points.length - 2], len, s, d)]))));
        points.push(new PPoint3d((i + .5) * len, len * k, 0, new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - 1], len, sc, dc), new FLinkage(points[points.length - 2], len, s, d)]))));
    }
    return points;
}

function rope(links: number, len = 1, s = 1, d = 1): PPoint3d[] {
    var points = [new PPoint3d(0, 0, 0, new noPhysics())];
    for (var i = 1; i < links; i++) {
        points.push(new PPoint3d(i * len, 0, 0, new GravityDecorator(0, g, 0, new PLinkage(points[points.length - 1], len, s, d))));
    }
    return points;
}

function stiffRopeNoAngle(links: number, len = 1, slen = 1, s = 1, d = 1, stiff = 1, stiffd = 1): PPoint3d[] {
    var points = [new PPoint3d(0, 0, 0, new noPhysics()), new PPoint3d(len, 0, 0, new noPhysics())];
    for (var i = 2; i < links; i++) {
        points.push(new PPoint3d(i * len, 0, 0, new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - 1], len, s, d), new FLinkage(points[points.length - 2], slen, stiff, stiffd)]))));
    }
    return points;
}
function stiffRope(links: number, len = 1, s = 1, d = 1, stiff = 1, stiffd = 1): PPoint3d[] {
    var points = [new PPoint3d(0, 0, 0, new noPhysics())];
    for (var i = 1; i < links; i++) {
        points.push(new PPoint3d(i * len, 0, 0, new GravityDecorator(0, g, 0, new PNeut(0, 0, 0, 1, [new FLinkage(points[points.length - 1], len, s, d), new FAngleLinkage(points[points.length - 1], defaultAngle, defaultAngle, stiff, stiffd)]))));
        if (i > 1) {
            ((points[points.length - 2].p as GravityDecorator).p as PNeut).f.push(new FAngleLinkage(points[points.length - 1], defaultAngle, defaultAngle, stiff, stiffd))
        }
    }
    return points;
}

export {
    Physics, noPhysics, Angle, QuaternionAngle,
    rope, PVel, PNeut, PLinkage, GravityDecorator, //CombinePhysics,
    TriBridge, Tube, FLinkage, stiffRope, Impulse, Hinge,
}
