#ifndef GDTOSVG_H
#define GDTOSVG_H
void svgImageCreate(FILE *svgout,int width, int height);

void svgEmptyRectangle(FILE *f, int x1, int y1, int x2, int  y2);
void svgFilledRectangle(FILE *f, int x1, int y1, int x2, int  y2,char *color);
void svgFilledRectangleWhiteBorder(FILE *f, int x1, int y1, int x2, int  y2,char *color);

void svgRectangle(FILE *f, int x1, int y1, int x2, int  y2,char *color);
void svgFilledRectangleNoBorder(FILE *svgout, int x1, int y1, int x2, int  y2,char *color);
void svgLine(FILE *f, int x1, int y1, int x2, int  y2,char *color);
void svgStringUp( FILE *f,int sizefont,int x1,int y1,char *st,char *color);
void svgString( FILE *f,int sizefont,int x1,int y1,char *st,char *color);
void svgDashedLine(FILE *svgout, int x1, int y1, int x2, int  y2,char *color);
void svgCircle(  FILE *svgout,int x1,int y1,int r, char *color);
void svgCircleJS(  FILE *svgout,int x1,int y1,int r, char *color,int nbp,double d,int g,char *name);
void svgCircleJS2(  FILE *svgout,int x1,int y1,int r, char *color,double d,double p);
void svgImageSetPixel(FILE *f, int x1, int y1, char *color) ;
void svgHorLine(FILE *svgout, int x1, int x2, int  y,char *color);
void svgVertLine(FILE *svgout, int x, int y1, int  y2,char *color);
void svgCircleBorder(  FILE *svgout,int x1,int y1,int r, char *color);
#endif