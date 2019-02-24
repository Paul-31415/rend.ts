import * as PIXI from "pixi.js";
import 'pixi-display';
import { traceOrtho, exampleScene, Ray, V3D, tracePersp, ColorImage, traceRefine } from "./render";
import { rope, PVel, PNeut, PLinkage, GravityDecorator, TriBridge, Tube, FLinkage, stiffRope, Angle, QuaternionAngle } from "./PixelArt";
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

//var rope1 = Tube(tubeC, 60, tubeL, 10000, 1, 10000, 1, 100);
var rope1 = stiffRope(3, tubeL, 1000, 1, 100000, 10000);

rope1[0].a = rope1[0].a.plusV([0, 0, 1], 0.01);
console.log((rope1[0].a as QuaternionAngle).axisAngle());
console.log((rope1[0].a as QuaternionAngle).dot(rope1[0].a as QuaternionAngle));
console.log(new QuaternionAngle(.5, .5, .5, .5).axisAngle());
console.log(new QuaternionAngle(1, 0, 1, 0).times(new QuaternionAngle(1, 0.5, .5, 0.75)));


//rope1[rope1.length - 1].p.vy = 1;

//var Rope1 = new PIXI.mesh.Rope(PIXI.Texture.fromImage("./static/line.png"), rope1);
//app.stage.addChild(Rope1);


for (var i = 1; i < rope1.length; i++) {
    //(rope1[i].p as PVel).vz = 0.1 * i * (i) / rope1.length * i / rope1.length;//-100;
    //(rope1[i].p as PVel).vx = 0.1 * i * (i) / rope1.length * i / rope1.length;//-100;
    //(rope1[i].p as PNeut).mass = 500 - i;
    //rope1[i].x = 0;//i % 2;
    //rope1[i].y = i;//i % 2;
    ((rope1[i].p as GravityDecorator).p as PNeut).I = 1000000;
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

const dscale = 0.0001;

var renderer = PIXI.autoDetectRenderer(800, 800);

var t = 0;

const ovc = 10;


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


function update(delta: number) {
    Rscene.clear();
    //rope1[0].y = mousePosition[1];

    Rscene.moveTo(0, 0);
    for (var j = 0; j < ovc; j++) {
        //((rope1[400].p as GravityDecorator).p as PLinkage).l += 0.001;

        t += delta / ovc;
        rope1[0].z = 0 * Math.sin(t / 14);

        /*for (var i = tubeC; i < rope1.length; i++) {
            for (var n = 0; n < ((rope1[i].p as GravityDecorator).p as PNeut).f.length; n++) {
                (((rope1[i].p as GravityDecorator).p as PNeut).f[n] as FLinkage).l *= 1 - delta / ovc * (i / rope1.length * i / rope1.length * i / rope1.length / 1000);
            }
        }*/
        //(rope1[499].p as GravityDecorator).gx += Math.random() - .5;
        //(rope1[499].p as GravityDecorator).gy += Math.random() - .5;//(rope1[499].p as PVel).vy;
        //point(t / 14, t / 51, t / 40)

        for (var i = 0; i < rope1.length; i++) {
            rope1[i].update(dscale * delta);

        }
        for (var i = rope1.length; i > 0; i--) {
            rope1[i - 1].update(dscale * delta);
        }
    }
    for (var i = 0; i < rope1.length; i++) {
        const z = Math.floor(rope1[i].z * 255 / 500 + 127);
        Rscene.lineStyle(2, z * 0x010101, 0.8);
        Rscene.lineTo(rope1[i].x, rope1[i].y);


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

}


async function startGame() {
    //    io.Socket.emit("testConsoleLog", "Hello from client");
    app.ticker.add(update);

}
startGame();
