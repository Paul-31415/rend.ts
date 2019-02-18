import * as express from "express";
import * as http from "http";
//import * as socket from "socket.io";
import * as path from "path";


const appRoot: string = path.join(__dirname, "../");

const app = express();
const httpServer = new http.Server(app);
//const io = socket(httpServer);

const port: number = 8000;

app.use("/static", express.static(__dirname + "/static"));
app.use("/", function(_req: express.Request, res: express.Response) {
    res.sendFile(__dirname + "/index.html");
});

httpServer.listen(port, function() {
    console.log("listening at port " + port);
});






/*
io.on('connection', function(socket) {
    socket.on("testConsoleLog", function(data) {
        console.log(data);
    });



});
*/
