import * as PIXI from "pixi.js";
import * as CTRL from "./controls";


import 'pixi-display';
import { traceOrtho, exampleScene, Ray, V3D, tracePersp, ColorImage, traceRefine } from "./render";
import { rope, PVel, PNeut, PLinkage, GravityDecorator, TriBridge, Tube, FLinkage, stiffRope, Angle, QuaternionAngle, FAngleLinkage, FStiffLinkage } from "./PixelArt";
import { test } from "mocha";
//import { io } from 'socket.io-client';


const app = new PIXI.Application({
    resolution: window.devicePixelRatio || 1,
    autoResize: true,
    width: window.innerWidth,
    height: window.innerHeight
});


document.body.appendChild(app.view);



const style = new PIXI.TextStyle({
    fill: "#ffffff"
});

const exampleText = new PIXI.Text("This is a test", style);
exampleText.position.x = 100;
exampleText.position.y = 100;
exampleText.anchor.x = 0.5;
exampleText.anchor.y = 0.5;


//app.stage.addChild(exampleText);

//var engine = new Engine();

const triangle = new PIXI.Graphics();
triangle.beginFill(0x20ff0000);
triangle.moveTo(-50, -50);
triangle.lineTo(50, -50);
triangle.lineTo(0, 36);
triangle.lineTo(-50, -50);
triangle.endFill();

triangle.position.x = 200;
triangle.position.y = 200;

//app.stage.addChild(triangle);

function mod(x: number, m: number): number {
    return ((x % m) + m) % m;
}
function round(x: number, v: number): number {
    return Math.floor(x / v) * v;
}

var ss = 2.5;

const Rsphere = new PIXI.Graphics();
for (var y = -50; y < 50; y += ss + ss + ss + ss) {

    Rsphere.moveTo(0, y);
    for (var x = -50; x < 50; x += ss + ss + ss + ss) {
        var color: number = 0;

        if (x * x + y * y < 40 * 40) {
            var z = (Math.sqrt(40 * 40 - x * x - y * y));
            color = round(Math.floor(mod(x / (z + 0.01), 1) * 256), 128) * 256 + round(mod(y / (z + 0.01), 1) * 256, 128);

        }

        Rsphere.lineStyle(ss, color, 1, 1);
        Rsphere.lineTo(x, y);

    }
}

Rsphere.position.x = 400;
Rsphere.position.y = 400;

//app.stage.addChild(Rsphere);

const samplesPerRound = 1;
const stepsPerRound = 50;
//var res: ColorImage = tracePersp(exampleScene, new Ray(new V3D(0, 0, -4), new V3D(-1.5, -1.5, 1)), new V3D(3, 0, 0), new V3D(0, 3, 0), 1 / 128 / 8, 1 / 128 / 8, 30, 1);

const tubeC = 6;
const tubeL = 10;

//var rope1 = Tube(tubeC, 60, tubeL, 10, 1, 1000, 1, 100);
var rope1 = stiffRope(24, [15, 0, 0], 20, 200, 10, 0, -100);

const r = rope1[0].a.plusV([1, 1, 1], 1);

//rope1[1].p.vy = 1;
//rope1[1].p.vz = 1;

//console.log((rope1[0].a as QuaternionAngle).axisAngle());
//console.log((rope1[0].a as QuaternionAngle).dot(rope1[0].a as QuaternionAngle));
//console.log(new QuaternionAngle(.5, .5, .5, .5).axisAngle()); 
//console.log(new QuaternionAngle(1,  0, 1, 0).times(new QuaternionAngle(1, 0.5, .5, 0.75)));

/*rope1[0].p = new PNeut(0, 0, 0, 1, [new FStiffLinkage(rope1[1], [-10, 0, 0], new QuaternionAngle(), (((rope1[1].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).stiffness
    , 10,
    (((rope1[1].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).dampening,
    (((rope1[1].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).dampening2, -100)], [0, 0, 0], 100);*/

//rope1 = rope1.reverse();
//rope1.pop();
//rope1 = rope1.reverse();

rope1[rope1.length - 1].p.vz = .1;


//rope1[0].a = new QuaternionAngle().plusV([0, 1, 1], 3);
//rope1[1].a = rope1[0].a.plusV([-1, -1, -1], 1);
//debugger;
//var Rope1 = new PIXI.mesh.Rope(PIXI.Texture.fromImage("./static/line.png"), rope1);
//app.stage.addChild(Rope1); 


for (var i = 1; i < rope1.length; i++) {
    //(rope1[i].p as PVel).vz = 0.1 * i * (i) / rope1.length * i / rope1.length;//-100;
    //(rope1[i].p as PVel).vx = 0.1 * i * (i) / rope1.length * i / rope1.length;//-100;
    //(rope1[i].p as PNeut).mass = 500 - i;
    //rope1[i].x = 0;//i % 2;
    //rope1[i].y = i;//i % 2;
    ((rope1[i].p as GravityDecorator).p as PNeut).I = 25 / 1;
    ((rope1[i].p as GravityDecorator).p as PNeut).mass = 1 / 1;
	/*if (i > 0) {
                (rope1[i].p as GravityDecorator).gy = Math.sin((i / 50) + t * t * 0.01);
                (rope1[i].p as GravityDecorator).gx = Math.cos((i / 50) + t * t * 0.01);
            }*/
}
//(rope1[499].p as PVel).vy = -30;
//(rope1[249].p as PVel).vy = 30;


var Rscene: PIXI.Graphics = new PIXI.Graphics();
Rscene.x = 400;
Rscene.y = 400;

Rscene.addChild(exampleText);

app.stage.addChild(Rscene);

const speed = 20;



const dscale = 0.1 * speed;

var renderer = PIXI.autoDetectRenderer(800, 800);

var t = -10;


const ovc = 2 * speed;


function point(rx = 0, ry = 0, rz = 0) {
    const r = tubeL / (2 * Math.sin(Math.PI / tubeC));
    for (var i = 0; i < tubeC; i++) {
        const a = r * Math.sin(i * Math.PI * 2 / tubeC + rx);
        const b = r * Math.cos(i * Math.PI * 2 / tubeC + rx);


        rope1[i].x = 0;
        rope1[i].y = a;
        rope1[i].z = b;

        const tmpx = rope1[i].x * Math.cos(ry) - rope1[i].z * Math.sin(ry);
        rope1[i].z = rope1[i].x * Math.sin(ry) + rope1[i].z * Math.cos(ry);
        rope1[i].x = tmpx;

        const tmpy = rope1[i].y * Math.cos(rz) - rope1[i].x * Math.sin(rz);
        rope1[i].x = rope1[i].y * Math.sin(rz) + rope1[i].x * Math.cos(rz);
        rope1[i].y = tmpy;

    }
}

const alen = 10;
var unsimulatedTime = 0;

var zflex = 0;
var yflex = 0;
const ctrls = new CTRL.Key4Control("KeyW",
    "KeyA",
    "KeyS",
    "KeyD");


function doSimStep(): number {
    //for (var j = 0; j < ovc; j++) {
    //((rope1[400].p as GravityDecorator).p as PLinkage).l += 0.001;
    const delta = 1;
    t += delta / ovc;
    /*
    if (t > 125 && t < 9000) {
        for (var i = 1; i < rope1.length; i++) {
            (((rope1[i].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).angle = (((rope1[i].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).angle.plusV([0, 0, (t % 100 > 50) ? .001 : -.001], dscale * delta / ovc) as QuaternionAngle;
        }
    }
    if (t > 10000) {
        for (var i = 1; i < rope1.length; i++) {
            //(((rope1[i].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).angle = (((rope1[i].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).angle.plusV([0, 0.001, 0], dscale * delta / ovc) as QuaternionAngle;
        }
    }//*/



    const v = ctrls.value();
    for (var i = 1; i < rope1.length; i++) {
        (((rope1[i].p as GravityDecorator).p as PNeut).f[0] as FStiffLinkage).angle = new QuaternionAngle().plusV([0, v[0], v[1]], 0.1) as QuaternionAngle;
    }





    rope1[0].a = new QuaternionAngle();
    rope1[0].p.va = [0, 0, 0];
    //rope1[0].p.va = [Math.cos(t / 160) / 10, 0, Math.sin(t / 160) / 10];
    //rope1[0].p.va[0] += Math.cos(t / 160) / 400;
    //rope1[0].p.va[1] += Math.sin(t / 160) / 400;

    rope1[0].x = 0;
    rope1[0].y = 0;
    rope1[0].z = 0;
    rope1[0].p.vx = 0;
    rope1[0].p.vy = 0;
    rope1[0].p.vz = 0;

    /*for (var i = tubeC; i < rope1.length; i++) {
        for (var n = 0; n < ((rope1[i].p as GravityDecorator).p as PNeut).f.length; n++) {
            (((rope1[i].p as GravityDecorator).p as PNeut).f[n] as FLinkage).l *= 1 - delta / ovc * (i / rope1.length * i / rope1.length * i / rope1.length / 1000);
        }
    }*/
    //(rope1[499].p as GravityDecorator).gx += Math.random() - .5;
    //(rope1[499].p as GravityDecorator).gy += Math.random() - .5;//(rope1[499].p as PVel).vy;
    //point(t / 14, t / 51, t / 40)

    for (var i = 0; i < rope1.length; i++) {
        rope1[i].update(dscale * delta / ovc);

    }
    for (var i = rope1.length; i > 0; i--) {
        rope1[i - 1].update(dscale * delta / ovc);
    }
    return delta / ovc;
}
function render(): void {
    Rscene.moveTo(rope1[0].x, rope1[0].y);
    for (var i = 0; i < rope1.length; i++) {
        const vs = 10;
        const z = Math.floor(rope1[i].z * 255 / 500 + 127);
        Rscene.lineStyle(2, z * 0x010101, 0.8);
        Rscene.lineTo(rope1[i].x, rope1[i].y);

        if (i > 0 && true) {
            const X = -200 + 20 * i;
            const Y = -200;
            Rscene.moveTo(X, Y);
            const r = rope1[i].a.minus(rope1[i - 1].a);

            Rscene.lineStyle(0.5, 0x0000ff, 1);
            const x = r.apply(alen, 0, 0);
            Rscene.lineTo(X + x[0], Y + x[1]);
            Rscene.moveTo(X, Y);


            Rscene.lineStyle(0.5, 0x00ff00, 1);
            const y = r.apply(0, alen, 0);
            Rscene.lineTo(X + y[0], Y + y[1]);
            Rscene.moveTo(X, Y);


            Rscene.lineStyle(0.5, 0xff0000, 1);
            const Z = r.apply(0, 0, alen);
            Rscene.lineTo(X + Z[0], Y + Z[1]);
            Rscene.moveTo(X, Y);


            const tq = (r as QuaternionAngle).axisAngle();
            Rscene.lineStyle(0.5, Math.floor(tq[2] * vs + 128) * 0x010100, 1);
            Rscene.lineTo(X + tq[0] * vs, Y + tq[1] * vs);
            Rscene.moveTo(X, Y + 100);

            Rscene.lineStyle(1, 0xffff00, 1);
            Rscene.lineTo(X, Y + 100 + tq[0] * vs);
            Rscene.moveTo(X + 3, Y + 100);
            Rscene.lineTo(X + 3, Y + 100 + tq[1] * vs);
            Rscene.moveTo(X + 6, Y + 100);
            Rscene.lineTo(X + 6, Y + 100 + tq[2] * vs);

            Rscene.lineStyle(1, 0x0000ff, 1);
            Rscene.moveTo(X, Y + 100 + vs * Math.PI);
            Rscene.lineTo(X - 10, Y + 100 + vs * Math.PI);
            Rscene.moveTo(X, Y + 100 - vs * Math.PI);
            Rscene.lineTo(X - 10, Y + 100 - vs * Math.PI);

            const s = (r as QuaternionAngle).dot(r as QuaternionAngle);
            Rscene.moveTo(X + 10, Y + 100);
            Rscene.lineStyle(1, 0xffffff, 1);
            Rscene.lineTo(X + 10, Y + 100 + s * vs);

            Rscene.moveTo(X, Y);



            const va1 = rope1[i].p.va;
            const va2 = rope1[i - 1].p.va;
            const va = [va2[0] - va1[0], va2[1] - va1[1], va2[2] - va1[2]];
            Rscene.lineStyle(0.5, Math.floor(va[2] * vs + 128) * 0x010001, 1);
            Rscene.lineTo(X + va[0] * vs, Y + va[1] * vs);


            Rscene.moveTo(rope1[i].x, rope1[i].y);
        }



        Rscene.lineStyle(0.5, z * 0x000001, 1);
        const x = rope1[i].a.apply(alen, 0, 0);
        Rscene.lineTo(rope1[i].x + x[0], rope1[i].y + x[1]);
        Rscene.moveTo(rope1[i].x, rope1[i].y);


        Rscene.lineStyle(0.5, z * 0x000100, 1);
        const y = rope1[i].a.apply(0, alen, 0);
        Rscene.lineTo(rope1[i].x + y[0], rope1[i].y + y[1]);
        Rscene.moveTo(rope1[i].x, rope1[i].y);


        Rscene.lineStyle(0.5, z * 0x010000, 1);
        const Z = rope1[i].a.apply(0, 0, alen);
        Rscene.lineTo(rope1[i].x + Z[0], rope1[i].y + Z[1]);
        Rscene.moveTo(rope1[i].x, rope1[i].y);

        Rscene.lineStyle(0.5, z * 0x010001, 1);
        const va = rope1[i].p.va;

        Rscene.lineTo(rope1[i].x + va[0] * vs, rope1[i].y + va[1] * vs);
        Rscene.moveTo(rope1[i].x, rope1[i].y);




        if (i >= tubeC) {
            for (var n = 0; n < ((rope1[i].p as GravityDecorator).p as PNeut).f.length; n++) {
                Rscene.lineStyle(1, 255 * (((n + 1) & 1) + 256 * (((n + 1) & 2) + 256 * ((n + 1) & 4))), 0.2);
                Rscene.lineTo((((rope1[i].p as GravityDecorator).p as PNeut).f[n] as FLinkage).link.x, (((rope1[i].p as GravityDecorator).p as PNeut).f[n] as FLinkage).link.y);

                Rscene.moveTo(rope1[i].x, rope1[i].y);

            }
        }
    }
    //exampleText.rotation += delta / 100 + delta * exampleText.rotation / 100;
    //Rsphere.rotation += delta / 100;
    //Rsphere.position.x = (Rsphere.position.x % 0.5) + 400.25;
    //Rsphere.position.y = (Rsphere.position.y % 0.5) + 400.125;

    //app.stage.removeChild(Rscene);
    //    traceRefine(exampleScene, res, stepsPerRound, samplesPerRound);

    //Rscene = res.get(.5);
    //app.stage.addChild(Rscene);
    exampleText.text = String(t);

}



function update(delta: number) {
    Rscene.clear();
    //rope1[0].y = mousePosition[1];


    Rscene.moveTo(0, 0);
    unsimulatedTime += delta;
    while (unsimulatedTime > 0) {
        unsimulatedTime -= doSimStep();
    }
    render();
}



async function startGame() {
    //    io.Socket.emit("testConsoleLog", "Hello from client");
    app.ticker.add(update);

}
startGame();
