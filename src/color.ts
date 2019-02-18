const bitsPerColor = 8;



abstract class Color {
    abstract getRGB(): number;
    abstract getR(): number;
    abstract getG(): number;
    abstract getB(): number;
    abstract getAlpha(): number;
    abstract scale(n: number): Color;
    abstract plus(c: Color): Color;
    abstract times(c: Color): Color;
}
function bound(n: number, l: number, h: number): number {
    return Math.max(l, Math.min(h, n));
}
function packColor(r: number, g: number, b: number): number {
    return (bound(Math.floor(r), 0, (1 << bitsPerColor) - 1) << (bitsPerColor * 2)) + (bound(Math.floor(g), 0, (1 << bitsPerColor) - 1) << bitsPerColor) + bound(Math.floor(b), 0, (1 << bitsPerColor) - 1);
}



class RGBColor extends Color {
    private static scalar = 255;

    constructor(public r: number,
        public g: number,
        public b: number) { super(); }


    public static convert(c: Color): RGBColor {
        return new RGBColor(c.getR(), c.getG(), c.getB());
    }


    public getRGB(): number {
        return packColor(this.r * RGBColor.scalar, this.g * RGBColor.scalar, this.b * RGBColor.scalar);
    }
    public getR(): number { return this.r; }
    public getG(): number { return this.g; }
    public getB(): number { return this.b; }
    public getAlpha(): number {
        return 1;
    }
    public scale(n: number): Color {
        return new RGBColor(this.r * n, this.g * n, this.b * n);
    }
    public plus(c: Color): Color {
        var cc = RGBColor.convert(c);
        return new RGBColor(this.r + cc.r, this.g + cc.g, this.b + cc.b);
    }
    public times(c: Color): Color {
        var cc = RGBColor.convert(c);
        return new RGBColor(this.r * cc.r, this.g * cc.g, this.b * cc.b);
    }
}



export {
    RGBColor, Color
}
