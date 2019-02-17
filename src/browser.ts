import * as PIXI from "pixi.js";
import 'pixi-display';

const app = new PIXI.Application({
    resolution: window.devicePixelRatio || 1,
    autoResize: true,
    width: window.innerWidth,
    height: window.innerHeight
});

document.body.appendChild(app.view);



const style = new PIXI.TextStyle({
    fill: "#ff3328"
});

const exampleText = new PIXI.Text("This is a test", style);
exampleText.position.x = 100;
exampleText.position.y = 100;
exampleText.anchor.x = 0.5;
exampleText.anchor.y = 0.5;


app.stage.addChild(exampleText);

//var engine = new Engine();


function update(delta: number) {
    exampleText.rotation += delta / 100;

}

function startGame() {
    app.ticker.add(update);
}


startGame();
