syms a b c d e
bending0 = (a + c - 2*b)^2;
bending1 = (cos(a) + cos(c) - 2*cos(b))^2 + (sin(a) + sin(c) - 2*sin(b))^2;
bending2 = 6 + 2*cos(a-c) - 4*cos(b-a) - 4*cos(b-c);
simplify(bending1 - bending2)



bending_b = 6 + 2*cos(a-c) - 4*cos(b-a) - 4*cos(b-c);
bending_c = 6 + 2*cos(b-d) - 4*cos(c-b) - 4*cos(d-c);
bending_d = 6 + 2*cos(e-c) - 4*cos(d-c) - 4*cos(e-d);
L = bending_b + bending_c + bending_d;

g_b = simplify(diff(L,b));
g_c = simplify(diff(L,c));
g_d = simplify(diff(L,d));

H_cc = simplify(diff(g_c,c))
H_cb = simplify(diff(g_c,b))
H_cd = simplify(diff(g_c,d))

% =========================================================================

ben0 = @(d1,d2) (d1 - d2).^2;
ben1 = @(d1,d2) 2*(1 - cos(d1 - d2));

d1 = linspace(-pi,pi);
d2 = linspace(-pi,pi);
figure
clf
hold on

val1 = ben1(d1,d2');
col1 = (val1-min(val1(:)))/(max(val1(:))-min(val1(:)));
col1 = cat(3,ones(size(col1)), col1, col1);
q1 = surf(d1,d2,val1,col1);
q1.EdgeAlpha = 0;
q1.FaceAlpha = 0.5;

val0 = ben0(d1,d2');
col0 = (val0-min(val0(:)))/(max(val0(:))-min(val0(:)));
col0 = cat(3,col0, col0, ones(size(col0)));
q0 = surf(d1,d2,val0,col0);
q0.EdgeAlpha = 0;
q0.FaceAlpha = 0.5;

hold off

% =========================================================================

mem0 = @(d) d.^2;
mem1 = @(d) 2*(1-cos(d));
d = linspace(-3*pi,3*pi);
figure
clf
hold on
val1 = mem1(d);
q1 = plot(d,val1,'r');
val0 = mem0(d);
q0 = plot(d,val0,'b');
hold off