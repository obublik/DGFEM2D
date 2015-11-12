function sterbina
l1 = 0.1;
l2 = 1;
d = 0.1;
d2 = 0.5;
x = [0 l1 l1 l1+l2 l1+l2 2*l1+l2 2*l1+l2 0];
y = [0 0 -d2 -d2 0 0 d d];
data = [x',y'];
save sterbina data;