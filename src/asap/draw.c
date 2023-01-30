 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "asap.h"
#include "gdtosvg.h"

#include "wrapio.h"

#define NBCOLORS 10

#define SIZEOFLEGEND 100
#define SIZEOFSQ 10
#define MARGE 50
int SIZEOFGRAPHIC;
// for values beetween 0 and 1
//valeur petite=petite valeur de v=rouge


unsigned char getColor(double val, double maxi, double mini)
{

	//int i;
	//double min;
	//double v;

	if (val> 0.2)
		return (255);
	else
		if (val> 0.1)
			return (210);
		else
			if (val> 0.05)
				return (170);
			else
				if (val> 0.01)
					return (130);
			else
				if (val> 0.005)
					return (90);
			else
				if (val> 0.001)
				return (50);
			else
				return (0);
}


void CreateHeadersvg(FILE *svgout)
{
	fprintf(svgout, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
	fprintf(svgout, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n") ;
	fprintf(svgout, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(svgout, "<html><head>\n");
	fprintf(svgout, "  <style>\n");
	fprintf(svgout, "    .caption{\n");
	fprintf(svgout, "	font-size: 14px;\n");
	fprintf(svgout, "	font-family: Georgia, serif;\n");
	fprintf(svgout, "    }\n");
	fprintf(svgout, "    .tooltip{\n");
	fprintf(svgout, "	font-size: 12px;\n");
	fprintf(svgout, "    }\n");
	fprintf(svgout, "    .tooltip_bg{\n");
	fprintf(svgout, "	fill: white;\n");
	fprintf(svgout, "	stroke: black;\n");
	fprintf(svgout, "	stroke-width: 1;\n");
	fprintf(svgout, "	opacity: 0.85;\n");
	fprintf(svgout, "   }\n");
	fprintf(svgout, "  </style>\n");

	fprintf(svgout, "  <script type=\"text/javascript\">\n");


	fprintf(svgout, "	function init(evt)\n");
	fprintf(svgout, "	{\n");
	fprintf(svgout, "	    if ( window.svgDocument == null )\n");
	fprintf(svgout, "	    {\n");
	fprintf(svgout, "		svgDocument = evt.target.ownerDocument;\n");
	fprintf(svgout, "	    }\n");

	fprintf(svgout, "	    tooltip = svgDocument.getElementById('tooltip');\n");
	fprintf(svgout, "	    tooltip_bg = svgDocument.getElementById('tooltip_bg');\n");

	fprintf(svgout, "	}\n");

	fprintf(svgout, "	function ShowTooltip(evt, mouseovertext)\n");
	fprintf(svgout, "	{\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"x\",evt.clientX+11);\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"y\",evt.clientY+27);\n");
	fprintf(svgout, "	    tooltip.firstChild.data = mouseovertext;\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"visibility\",\"visible\");\n");

	fprintf(svgout, "	    length = tooltip.getComputedTextLength();\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"width\",length+8);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"x\",evt.clientX+8);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"y\",evt.clientY+14);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"visibility\",\"visibile\");\n");
	fprintf(svgout, "	}\n");

	fprintf(svgout, "	function HideTooltip(evt)\n");
	fprintf(svgout, "	{\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"visibility\",\"hidden\");\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"visibility\",\"hidden\");\n");
	fprintf(svgout, "	}\n");


	fprintf(svgout, "  </script>\n");
	fprintf(svgout, "</head>\n");
	fprintf(svgout, "<svg xmlns=\"http://www.w3.org/2000/svg\" onload=\"init(evt)\" ");
	fprintf(svgout, "width=\"%d\" height=\"%d\" >\n", SIZEOFGRAPHIC, SIZEOFGRAPHIC);

//<div style="width:320px; height:480px;"> OK FOR SCALING
//    <svg id="svg1" width="100%" height="100%" viewBox="0 0 640 480" preserveAspectRatio="xMaxYMax"></svg>
//</div>

	fprintf(svgout, "<g>\n"); // g for grouping all the graphic otherwise each square is an object
}



/*****************************/
void draw_heat_svg(FILE *f, struct DistanceMatrix mat, Node *zenode   , double min, double max)
{
	double ech, v;

	int i, j, ii,jj, c, aa, bb;
//char  *LavaColors[NBCOLORS]={"#FFFF00","#FFD700","#FFAF00","#FF8700","#F00000","#780000","#320000"};

	SIZEOFGRAPHIC = mat.n * SIZEOFSQ + (4 * MARGE);
	ech = SIZEOFSQ;


	CreateHeadersvg	(f);

	for (j = 0; j <	mat.n; j++)
		{
		ii = zenode[j].first_to_draw;
				fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"10\" onmousemove=\"ShowTooltip(evt, \'Seq%d:%s\')\" onmouseout=\"HideTooltip(evt)\"  style=\"writing-mode: tb;\"> S%d </text>\n",
				        (int)(ii*ech)+MARGE, 1, ii, mat.names[ii], j);

		}
	for (i = 0; i < mat.n; i++)
		{
		ii = zenode[i].first_to_draw;
		fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"10\"  onmousemove=\"ShowTooltip(evt, \'Seq%d:%s\')\" onmouseout=\"HideTooltip(evt)\"> S%d </text>\n",
				1, MARGE+(int)(ii*ech), ii, mat.names[ii], ii);
		for (j = 0; j < mat.n; j++)
			{
				jj=zenode[j].first_to_draw;
//				bb=ech + (jj*ech);
				v = mat.dist[ii][jj];
				c = getColor(v, max, min);

				aa =ii * ech;
				bb = jj * ech;
				fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill= \"rgb(%d,%d,%d)\" onmousemove=\"ShowTooltip(evt, \'Group %d d=%f\')\" onmouseout=\"HideTooltip(evt)\" />\n",
				        MARGE + aa, MARGE + bb, (int)ech, (int)ech, c, c,c,i + 1, mat.dist[ii][jj]);
			}
		}


	fprintf(f, "</g><rect class=\"tooltip_bg\" id=\"tooltip_bg\"\n");
	fprintf(f, "x=\"0\" y=\"0\" rx=\"4\" ry=\"4\"\n");
	fprintf(f, "width=\"55\" height=\"17\" visibility=\"hidden\"/>\n");
	fprintf(f, "<text class=\"tooltip\" id=\"tooltip\"\n");
	fprintf(f, "x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>\n");

	fprintf(f, "</svg></HTML>\n");

}
