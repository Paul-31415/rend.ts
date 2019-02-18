import * as PIXI from "pixi.js";
import 'pixi-display';
import { traceOrtho, exampleScene, Ray, V3D, tracePersp } from "./render";
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


app.stage.addChild(exampleText);

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

const Rscene = tracePersp(exampleScene, new Ray(new V3D(0, 0, -4), new V3D(-1.5, -1.5, 1)), new V3D(3, 0, 0), new V3D(0, 3, 0), 1 / 512, 1 / 512, 10, 100);




Rscene.position.x = 50;
Rscene.position.y = 50;


app.stage.addChild(Rscene);




function update(delta: number) {
    exampleText.rotation += delta / 100 + delta * exampleText.rotation / 100;
    Rsphere.rotation += delta / 100;
    Rsphere.position.x = (Rsphere.position.x % 0.5) + 400.25;
    Rsphere.position.y = (Rsphere.position.y % 0.5) + 400.125;
}


async function startGame() {
    //    io.Socket.emit("testConsoleLog", "Hello from client");
    app.ticker.add(update);

}
startGame();
