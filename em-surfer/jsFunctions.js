var bgR = 0;
var bgG = 0;
var bgB = 0;
var sX = 0;
var sY = 0;
var vginit = 0;
var vgcavity = "";
var vgprotrusion = "";
var vgflat = "";
var ceArray = new Array(25);

var newwindow;

function poptastic(url)
{
	newwindow=window.open(url,'name','height=400,width=200');
		if (window.focus) {newwindow.focus()}
}

function toHex(i)
{
	i = parseInt(i).toString(16);
	return i.length < 2 ? "0" + i : i;
}

function changeBG(event)
{
        var x = event.offsetX?(event.offsetX):event.pageX-document.getElementById("colorbar").offsetLeft;
        var y = event.offsetY?(event.offsetY):event.pageY-document.getElementById("colorbar").offsetTop;
	var c = getRGB(y * 500 + x - 1, 500, 30).toString();
	var temp = new Date();
	var days = 30;
	temp.setTime(temp.getTime() + days * 24 * 60 * 60 * 1000);
//        document.bgColor = '#' + getRGB(y * 350 + x - 1, 350, 21).toString();
//        document.bgColor = '#' + getRGB(y * 500 + x - 1, 500, 30).toString();
        document.bgColor = '#' + c;
	document.jmol.script('color background [x' + c + '];');
//	document.jmol.script('color background [x' + getRGB(y * 500 + x - 1, 500, 30).toString() + '];');
	document.body.style.color = (255 - (bgR * 0.299 + bgG * 0.587 + bgB * 0.114)) < 105 ? '#000000' : '#FFFFFF';
//        document.jmol.boxbgcolor.value = getRGB(y * 350 + x - 1, 350, 21);
//      alert(document.bgColor);
//	document.cookie = "bgcolor=" + c + ";expires=" + temp.toGMTString();
	document.cookie = "bgcolor=" + c + ";expires=" + temp.toGMTString();
	document.cookie = "txtcolor=" + ((255 - (bgR * 0.299 + bgG * 0.587 + bgB * 0.114)) < 105 ? '#000000' : '#FFFFFF') + ";expires=" + temp.toGMTString();
//	alert (document.body.style.color);
}

function displayCE(i)
{
//	ceWindow = window.open("", "cewindow", "height=600, width=800, resizeable=1, status=1, scrollbars=1");
	ceWindow = window.open("", "cewindow", "menubar=1, resizeable=1, status=1, scrollbars=1");
	ceWindow.focus();
	ceWindow.document.write("<html><head></head><body><pre>");
	ceWindow.document.write(ceArray[i]);
	ceWindow.document.write("</pre><form method=\"post\"><input type=\"button\" onclick=\"window.close()\" value=\"Close Window\" /></form>");
	ceWindow.document.write("</body></html>");
	ceWindow.document.close();
}

function init(c, t)
{
	document.bgColor = '#' + c;
	document.jmol.script('color background [x' + c + '];');
	document.body.style.color = t;
//	if (bgR != 0 && bgG != 0 && bgR != 0)
//	{
//		document.body.style.color = (255 - (bgR * 0.299 + bgG * 0.587 + bgB * 0.114)) < 105 ? '#000000' : '#FFFFFF';
//	}
}

function getRGB(x, width, height)
{
        var wid = width / 7;
        var col = x % width;
        var row = Math.floor(x / (wid * 7));
        var c = x % wid;
        var r = 0;
        var g = 0;
        var b = 0;
        var h = 0;
        var pct = (255 / wid) * c;

        if (col >= wid * 6)
        {
                r = 255 - pct;
                g = 255 - pct;
                b = 255 - pct;
        }
        else
        {
                h = 255 - pct;
                r = col < wid ? 255 : col < wid * 2 ? h : col < wid * 4 ? 0 : col < wid * 5 ? pct : 255;
                g = col < wid ? pct : col < wid * 3 ? 255 : col < wid * 4 ? h : 0;
                b = col < wid * 2 ? 0 : col < wid * 3 ? pct : col < wid * 5 ? 255 : h;

                if (row < (height / 2))
                {
                        var base = 255 - (255 * 2 / height) * row;
                        r = base + (r * row * 2 / height);
                        g = base + (g * row * 2 / height);
                        b = base + (b * row * 2 / height);
                }
                else if (row > (height / 2))
                {
                        var base = (height - row) / (height / 2);
                        r = r * base;
                        g = g * base;
                        b = b * base;
                }
        }
	bgR = r;
	bgG = g;
	bgB = b;
        return toHex(r) + toHex(g) + toHex(b);
}

function spinX()
{
	if (sX == 0)
	{
		document.jmol.script('spin off; set spin Y 0; set spin X 0; set spin Y 50; spin on;');
		sY = 0;
		sX = 1;
	}
	else
	{
		document.jmol.script('spin off');
		sX = 0;
	}
}

function spinY()
{
	if (sY == 0)
	{
		document.jmol.script('spin off; set spin Y 0; set spin X 0; set spin X 50; spin on;');
		sX = 0;
		sY = 1;
	}
	else
	{
		document.jmol.script('spin off');
		sY = 0;
	}
}

function spin360()
{
	document.jmol.script('spin off; move 0 360 0 0 0 0 0 0 4; move 360 0 0 0 0 0 0 0 4;');
}
