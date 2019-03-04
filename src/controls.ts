//import { Engine, Input, GamepadConnectEvent, GameEvent, GamepadButtonEvent, GamepadAxisEvent, } from 'excalibur';

//const engine = new Engine();
/*
engine.input.gamepads.on('connect', (cee: GameEvent<any> | undefined) => {
    if (cee instanceof GamepadConnectEvent) {
        const ce = cee as GamepadConnectEvent;



        console.log('Gamepad connected', ce);
        ce.gamepad.on('button', (ev: GameEvent<any> | undefined) => {
            if (ev instanceof GamepadButtonEvent) {
                const be = ev as GamepadButtonEvent;
                if (be.button === Input.Buttons.Face1) {
                    console.log("facE!1");
                }
            }
        });

        ce.gamepad.on('axis', (ev: GameEvent<any> | undefined) => {
            if (ev instanceof GamepadAxisEvent) {
                const ae = ev as GamepadAxisEvent;
                if (ae.axis === Input.Axes.LeftStickX && ae.value > 0.5) {
                    console.log("rightytrews");
                }
            }

        });
    }
});
*/



interface BooleanControl {
    value(): boolean;
}
interface Number1DControl {
    value(): number;
}
interface Number2DControl {
    value(): number[];
}


class Key4Control implements Number2DControl {

    public keyUp: KeyControl;
    public keyLeft: KeyControl;
    public keyDown: KeyControl;
    public keyRight: KeyControl;



    constructor(keyUp: string,
        keyLeft: string,
        keyDown: string,
        keyRight: string
    ) {
        this.keyUp = new KeyControl(keyUp);
        this.keyLeft = new KeyControl(keyLeft);
        this.keyDown = new KeyControl(keyDown);
        this.keyRight = new KeyControl(keyRight);
    }
    value(): number[] {
        var v = [0, 0];
        if (this.keyUp.value()) {
            v[0] += 1;
        }
        if (this.keyDown.value()) {
            v[0] -= 1;
        }
        if (this.keyLeft.value()) {
            v[1] += 1;
        }
        if (this.keyRight.value()) {
            v[1] -= 1;
        }
        return v;



    }
}



class KeyControl implements BooleanControl {
    private val = false;
    constructor(public key: string) {
        document.addEventListener("keydown", (ke: KeyboardEvent) => {
            if (ke.code == this.key) {
                this.val = true;
            }
        });
        document.addEventListener("keyup", (ke: KeyboardEvent) => {
            if (ke.code == this.key) {
                this.val = false;
            }
        });

    }
    value(): boolean {

        return this.val;
    }

}



export {

    BooleanControl, KeyControl,
    Number1DControl,
    Number2DControl, Key4Control
}
