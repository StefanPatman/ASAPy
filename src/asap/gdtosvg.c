/*translate every gd fns ton svg using quite same syntax*/
#include <stdio.h>
#include "wrapio.h"

void svgImageCreate(FILE *svgout,int width, int height)
{

//fprintf(svgout,"<svg xmlns=\"http://www.w3.org/2000/svg\" onload=\"init(evt)\" ");
fprintf(svgout,"<svg xmlns=\"http://www.w3.org/2000/svg\"  ");
fprintf(svgout,"width=\"%d\" height=\"%d\" >\n", width, height);
}


void svgFilledRectangle(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{

 fprintf(svgout, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",x1,y1,x2,y2);
 fprintf(svgout, " style=\"fill:%s;stroke:black;\" />\n",color);
}

void svgFilledRectangleWhiteBorder(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{

 fprintf(svgout, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",x1,y1,x2,y2);
 fprintf(svgout, " style=\"fill:%s;stroke:white;\" />\n",color);
}
void svgEmptyRectangle(FILE *svgout, int x1, int y1, int x2, int  y2)
{

 fprintf(svgout, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",x1,y1,x2,y2);
 fprintf(svgout, " style=\"stroke:black;fill:none\" />\n");
}
void svgFilledRectangleNoBorder(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{

 fprintf(svgout, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",x1,y1,x2,y2);
 fprintf(svgout, " style=\"fill:%s;stroke:%s;\" />\n",color,color);
}

void svgRectangle(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{
fprintf(svgout, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",x1,y1,x2,y2);
fprintf(svgout, " style=\"stroke:%s;\" />\n",color);
}

void svgLine(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{
 fprintf(svgout, " <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"stroke:%s;stroke-width:1\" />\n",x1,y1,x2,y2,color);
}

void svgHorLine(FILE *svgout, int x1, int x2, int  y,char *color)
{
 fprintf(svgout, " <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"stroke:%s;stroke-width:1\" />\n",x1,y,x2,y,color);
}

void svgVertLine(FILE *svgout, int x, int y1, int  y2,char *color)
{
 fprintf(svgout, " <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"stroke:%s;stroke-width:1\" />\n",x,y1,x,y2,color);
}


void svgStringUp( FILE *svgout,int sizefont,int x1,int y1,char *st,char *color)
{
 fprintf(svgout, "  <text x=\"%d\" y=\"%d\" style=\"writing-mode: tb;color:%s;font-size  : %d;\" >%s</text>\n",x1,y1,color,sizefont,st);
}


void svgString( FILE *svgout,int sizefont,int x1,int y1,char *st,char *color)
{
 fprintf(svgout, "  <text x=\"%d\" y=\"%d\" style=\"color:%s;font-size :%d;font-family = monospace;\" >%s</text>\n",x1,y1,color,sizefont,st);

}
void svgCircle(  FILE *svgout,int x1,int y1,int r, char *color)
	{
  fprintf(svgout, " <circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" />\n",x1,y1,r,color);
	}
void svgCircleBorder(  FILE *svgout,int x1,int y1,int r, char *color)
	{
  fprintf(svgout, " <circle cx=\"%d\" cy=\"%d\" r=\"%d\" style=\"fill:%s;stroke:black;\" />\n",x1,y1,r,color);
	}

void svgCircleTitle(  FILE *svgout,int x1,int y1,int r, char *color,char *leg)
	{
  fprintf(svgout, " <circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" title=\"%s\" />\n",x1,y1,r,color,leg);
	}

void svgCircleJS(  FILE *svgout,int x1,int y1,int r, char *color,int nbp,double d,int g,char *objname)
	{
  fprintf(svgout, " <circle id=\"%s%d\" onclick=\"circle_click(evt,%d,%f,%d)\" cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" />\n",objname,nbp,nbp,d,g,x1,y1,r,color);
	}

void svgCircleJS2(  FILE *svgout,int x1,int y1,int r, char *color,double d,double p)
	{
  fprintf(svgout, " <circle onclick=\"circle_proba(evt,%f,%f)\" cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" />\n",d,p,x1,y1,r,color);
	}

void svgImageSetPixel(FILE *svgout, int x1, int y1, char *color)
{
		svgFilledRectangle(svgout, x1, y1, 1, 1,color)	;
}

void svgDashedLine(FILE *svgout, int x1, int y1, int x2, int  y2,char *color)
{
 fprintf(svgout, "    <line       x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-dasharray=\"5, 10\"/>\n",x1,y1,x2,y2,color);



 }
