boolean three_dimentions = false;
// gravitational constant for unit of mass
double gravity = 0.01;
long rand_seed = 1;
double win_border = 0.1; // border of screen to ensure is empty, relative to smallest window edge
// TODO: set to 0 if not random - same with ball_min
double max_size      = 40;  // max factor to zoom out by, can be set to 0 for infinite - only applies to random generation

// ball measurements are diameter
int    ball_count    = 3;
double ball_min      = 0.1; // relative to smallest window edge
double ball_max_mult = 2.5; // max ball size is this*min
double ball_speed    = 0.1; // max inital ball speed in relative-to-smallest-window-edge/frame

long cur_rand_seed = rand_seed;
double[][] balls = new double[ball_count][11]; // size, x, y, z, vx, vy, vz, done for this frame, last substep x, last substep y, last substep z
ArrayList<Double>[]  ballXs = new ArrayList[ball_count];
ArrayList<Double>[]  ballYs = new ArrayList[ball_count];
ArrayList<Double>[]  ballZs = new ArrayList[ball_count];
ArrayList<Double>[] ballVXs = new ArrayList[ball_count];
ArrayList<Double>[] ballVYs = new ArrayList[ball_count];
ArrayList<Double>[] ballVZs = new ArrayList[ball_count];
ArrayList<Double>[] ballAXs = new ArrayList[ball_count];
ArrayList<Double>[] ballAYs = new ArrayList[ball_count];
ArrayList<Double>[] ballAZs = new ArrayList[ball_count];

/*{{{ double rand()*/
// https://github.com/Qqwy/SimpleRNG
long max_32 = 1L << 32;
double rand() {
	cur_rand_seed ^= (cur_rand_seed << 13)%max_32;
	cur_rand_seed ^= (cur_rand_seed >> 17)%max_32;
	cur_rand_seed ^= (cur_rand_seed <<  5)%max_32;
	return (double)cur_rand_seed;
}
/*}}}*/

/*{{{ double interp(ArrayList<Double> list, double x)*/
double interp(ArrayList<Double> list, double x) {
	if(x == 0) return list.get(1);
	if(x/list.get(0) >= 1) return list.get(list.size()-1);
	double index = (x/list.get(0))*(list.size()-2)+1;
	return (1-index%1)*list.get((int)index)+(index%1)*list.get((int)index+1);
}
/*}}}*/

/*{{{ PVector rotate_v(float x, float y, float z)*/
int last_x, last_y;
float r[] = new float[3];
// rotate about x then y then z
PVector rotate_v(float x, float y, float z) {
	if(!three_dimentions) return new PVector(x, y, z);
	y = -y;
	PVector v;
	v = (new PVector(z,   y)).rotate(r[0]); // y and z
	y = v.y;
	v = (new PVector(x, v.x)).rotate(r[1]); // z and x
	z = v.y;
	v = (new PVector(y, v.x)).rotate(r[2]); // y and x
	return new PVector(v.y, -v.x, z);
}
PVector rotate_v(PVector v) { return rotate_v(v.x, v.y, v.z); }
float rotate_x(float x, float y, float z) { return rotate_v(x, y, z).x; }
float rotate_y(float x, float y, float z) { return rotate_v(x, y, z).y; }
float rotate_z(float x, float y, float z) { return rotate_v(x, y, z).z; }
float rotate_x(PVector v) { return rotate_v(v.x, v.y, v.z).x; }
float rotate_y(PVector v) { return rotate_v(v.x, v.y, v.z).y; }
float rotate_z(PVector v) { return rotate_v(v.x, v.y, v.z).z; }
/*}}}*/

/*{{{ void setup()*/
void setup() {
	for(int i = 0; i < ball_count; i++) {
		ballXs[i]  = new ArrayList<Double>();
		ballYs[i]  = new ArrayList<Double>();
		ballZs[i]  = new ArrayList<Double>();
		ballVXs[i] = new ArrayList<Double>();
		ballVYs[i] = new ArrayList<Double>();
		ballVZs[i] = new ArrayList<Double>();
		ballAXs[i] = new ArrayList<Double>();
		ballAYs[i] = new ArrayList<Double>();
		ballAZs[i] = new ArrayList<Double>();
		// size
		balls[i][0] = ball_min+((ball_max_mult-1)*ball_min)*(rand()/max_32);
		// position (edge of circle at -1 to 1)
		balls[i][1] = (1-balls[i][0]/2)*(2*(rand()/max_32)-1);
		balls[i][2] = (1-balls[i][0]/2)*(2*(rand()/max_32)-1);
		balls[i][3] = (1-balls[i][0]/2)*(2*(rand()/max_32)-1);
		// velocity
		balls[i][4] = ball_speed*(rand()/max_32);
		balls[i][5] = balls[i][4]/3;
		balls[i][6] = balls[i][4]/3;
		balls[i][4] = balls[i][4]/3;
		if(!three_dimentions) {
			balls[i][3] = 0;
			balls[i][6] = 0;
			balls[i][4] *= 3/2.0;
			balls[i][5] *= 3/2.0;
		}
	}
	size(600, 600, P2D);
}
/*}}}*/

void keyPressed() {
	if(key == ' ') {
		rand_seed++; cur_rand_seed = rand_seed;
		setup();
	}
}

void draw() {

/*{{{ camera rotation*/
	if(three_dimentions && mousePressed) {
		PVector x = rotate_v(1, 0, 0);
		PVector y = rotate_v(0, 1, 0);
		PVector z = rotate_v(0, 0, 1);
		PVector v;
		float inc;
		inc = 3.0/2*PI*(mouseX-last_x)/min(width, height);
		v = (new PVector(x.x, x.z)).rotate(inc);
		x.x = v.x; x.z = v.y;
		v = (new PVector(y.x, y.z)).rotate(inc);
		y.x = v.x; y.z = v.y;
		v = (new PVector(z.x, z.z)).rotate(inc);
		z.x = v.x; z.z = v.y;
		inc = -3.0/2*PI*(mouseY-last_y)/min(width, height);
		v = (new PVector(x.z, x.y)).rotate(inc);
		x.z = v.x; x.y = v.y;
		v = (new PVector(y.z, y.y)).rotate(inc);
		y.z = v.x; y.y = v.y;
		v = (new PVector(z.z, z.y)).rotate(inc);
		z.z = v.x; z.y = v.y;
		// ensure x axis has no y component
		if(x.y == 0) r[2] = 0;
		else {
			// get angle to move from correct x (no y component) to actual x
			r[2] = atan2(x.y, x.x);
			x = new PVector(sqrt(1-pow(x.z, 2)), 0, x.z);
			v = (new PVector(y.x, y.y)).rotate(-r[2]);
			y.x = v.x; y.y = v.y;
		}
		// ensure x axis has no z component
		// and that x axis is pointing right, only x rotation after (before?) this
		if(x.z == 0) {
			if(x.x > 0) r[1] = 0;
			else {
				r[1] = PI;
				x = new PVector(1, 0, 0);
				y = new PVector(-y.x, y.y, -y.z);
			}
		}
		else {
			r[1] = atan2(x.z, x.x);
			x = new PVector(1, 0, 0);
			v = (new PVector(y.x, y.z)).rotate(-r[1]);
			y.x = v.x; y.z = v.y;
		}
		// ensure y axis has no z component
		if(y.z == 0) r[0] = 0;
		// PI/2- is necessary because I did these calculations for y being up but it is down
		// y input and output are also negated in rotate_v
		else r[0] = PI/2-atan2(y.y, y.z);
	}
	last_x = mouseX;
	last_y = mouseY;
/*}}}*/

	background(255);
	translate(width/2, height/2);
	strokeWeight(0);
	fill(0);

/*{{{ camera center calculation*/
	// mean
	double mx = 0, my = 0, mz = 0;
	for(int i = 0; i < ball_count; i++) {
		mx += balls[i][1];
		my += balls[i][2];
		mz += balls[i][3];
	}
	mx /= ball_count; my /= ball_count; mz /= ball_count;
	for(int i = 0; i < ball_count; i++) {
		balls[i][1] -= mx;
		balls[i][2] -= my;
		balls[i][3] -= mz;
	}
	// least-squared-error mean, graph x=(some array) and y=sum_0^len of (array[i] - x)^2
	mx = 0; my = 0;
	double left_error, right_error;
	double min_x = Double.MAX_VALUE, max_x = -Double.MAX_VALUE;
	double min_y = Double.MAX_VALUE, max_y = -Double.MAX_VALUE;
	for(int i = 0; i < ball_count; i++) {
		PVector ball = rotate_v((float)balls[i][1], (float)balls[i][2], (float)balls[i][3]);
		min_x = Math.min(min_x, ball.x); max_x = Math.max(max_x, ball.x);
		min_y = Math.min(min_y, ball.y); max_y = Math.max(max_y, ball.y);
	}
	left_error = 0; right_error = 0;
	for(int i = 0; i < ball_count; i++) {
		float ball_x = rotate_x((float)balls[i][1], (float)balls[i][2], (float)balls[i][3]);
		left_error  += Math.pow(ball_x-min_x, 2);
		right_error += Math.pow(ball_x-max_x, 2);
	}
	mx = min_x+((right_error-left_error)/ball_count-Math.pow(max_x-min_x, 2))/(-2*(max_x-min_x));
	left_error = 0; right_error = 0;
	for(int i = 0; i < ball_count; i++) {
		float ball_y = rotate_y((float)balls[i][1], (float)balls[i][2], (float)balls[i][3]);
		left_error  += Math.pow(ball_y-min_y, 2);
		right_error += Math.pow(ball_y-max_y, 2);
	}
	my = min_y+((right_error-left_error)/ball_count-Math.pow(max_y-min_y, 2))/(-2*(max_y-min_y));
/*}}}*/

/*{{{ zoom*/
	double factor = 1;
	for(int i = 0; i < ball_count; i++) {
		PVector ball = rotate_v((float)balls[i][1], (float)balls[i][2], (float)balls[i][3]);
		factor = Math.max(factor, Math.max(Math.abs(ball.x-mx), Math.abs(ball.y-my))+balls[i][0]/2);
	}
	if(max_size != 0 && factor > max_size) {
		rand_seed++; cur_rand_seed = rand_seed;
		setup();
	}
	scale((float)(Math.min(width, height)/2/(factor/(1-2*win_border))));
/*}}}*/

	for(int i = 0; i < ball_count; i++) balls[i][7] = 0;
	boolean narrowed = false;
	for(int mult_step = 1; !narrowed && mult_step > 0; mult_step = mult_step << 1) {
		for(int i = 0; i < ball_count; i++) {
			if(balls[i][7] != 0) continue;
			 ballXs[i].clear();  ballXs[i].add(0d);  ballXs[i].add(balls[i][1]);
			 ballYs[i].clear();  ballYs[i].add(0d);  ballYs[i].add(balls[i][2]);
			 ballZs[i].clear();  ballZs[i].add(0d);  ballZs[i].add(balls[i][3]);
			ballVXs[i].clear(); ballVXs[i].add(0d); ballVXs[i].add(balls[i][4]);
			ballVYs[i].clear(); ballVYs[i].add(0d); ballVYs[i].add(balls[i][5]);
			ballVZs[i].clear(); ballVZs[i].add(0d); ballVZs[i].add(balls[i][6]);
			ballAXs[i].clear(); ballAXs[i].add(0d);
			ballAYs[i].clear(); ballAYs[i].add(0d);
			ballAZs[i].clear(); ballAZs[i].add(0d);
		}
		for(double i = 0; i < 1; i += 1f/mult_step) {
			for(int j = 0; j < ball_count; j++) {
				if(balls[j][7] != 0) continue;
				double ax = 0, ay = 0, az = 0;
				for(int k = 0; k < ball_count; k++) {
					if(k == j) continue;
					// max to prevent infinite acceleration when no collision checking
					double a = -gravity*balls[k][0]/Math.max(Math.pow(ball_min, 2),
						Math.pow(interp(ballXs[j], i)-interp(ballXs[k], i), 2)+
						Math.pow(interp(ballYs[j], i)-interp(ballYs[k], i), 2)+
						Math.pow(interp(ballZs[j], i)-interp(ballZs[k], i), 2)
					);
					double dx = interp(ballXs[j], i)-interp(ballXs[k], i);
					double dy = interp(ballYs[j], i)-interp(ballYs[k], i);
					double dz = interp(ballZs[j], i)-interp(ballZs[k], i);
					double d = Math.sqrt(dx*dx+dy*dy+dz*dz);
					ax += a*dx/d; ay += a*dy/d; az += a*dz/d;
				}
				ballXs[j].add(interp(ballXs[j], i)+(interp(ballVXs[j], i)+ax/mult_step/2)/mult_step); ballXs[j].set(0, i+1f/mult_step);
				ballYs[j].add(interp(ballYs[j], i)+(interp(ballVYs[j], i)+ay/mult_step/2)/mult_step); ballYs[j].set(0, i+1f/mult_step);
				ballZs[j].add(interp(ballZs[j], i)+(interp(ballVZs[j], i)+az/mult_step/2)/mult_step); ballZs[j].set(0, i+1f/mult_step);
				ballAXs[j].add(ax); ballAXs[j].set(0, i);
				ballAYs[j].add(ay); ballAYs[j].set(0, i);
				ballAZs[j].add(az); ballAZs[j].set(0, i);
			}
			for(int j = 0; j < ball_count; j++) {
				if(balls[j][7] != 0) continue;
				double ax = 0, ay = 0, az = 0;
				for(int k = 0; k < ball_count; k++) {
					if(k == j) continue;
					// max to prevent infinite acceleration when no collision
					double a = -gravity*balls[k][0]/Math.max(Math.pow(ball_min, 2),
						Math.pow(interp(ballXs[j], i+1f/mult_step)-interp(ballXs[k], i+1f/mult_step), 2)+
						Math.pow(interp(ballYs[j], i+1f/mult_step)-interp(ballYs[k], i+1f/mult_step), 2)+
						Math.pow(interp(ballZs[j], i+1f/mult_step)-interp(ballZs[k], i+1f/mult_step), 2)
					);
					double dx = interp(ballXs[j], i+1f/mult_step)-interp(ballXs[k], i+1f/mult_step);
					double dy = interp(ballYs[j], i+1f/mult_step)-interp(ballYs[k], i+1f/mult_step);
					double dz = interp(ballZs[j], i+1f/mult_step)-interp(ballZs[k], i+1f/mult_step);
					double d = Math.sqrt(dx*dx+dy*dy+dz*dz);
					ax += a*dx/d; ay += a*dy/d; az += a*dz/d;
				}
				ballVXs[j].add(interp(ballVXs[j], i)+(interp(ballAXs[j], i)+ax)/mult_step/2); ballVXs[j].set(0, i+1f/mult_step);
				ballVYs[j].add(interp(ballVYs[j], i)+(interp(ballAYs[j], i)+ay)/mult_step/2); ballVYs[j].set(0, i+1f/mult_step);
				ballVZs[j].add(interp(ballVZs[j], i)+(interp(ballAZs[j], i)+az)/mult_step/2); ballVZs[j].set(0, i+1f/mult_step);
			}
		}
		narrowed = true;
		for(int i = 0; i < ball_count; i++) {
			if(balls[i][7] != 0) continue;
			if(mult_step != 1 && (mult_step >= 64 ||
				Math.abs(1-(ballXs[i].get(ballXs[i].size()-1)-balls[i][1])/(balls[i][8 ]-balls[i][1]))+
				Math.abs(1-(ballYs[i].get(ballYs[i].size()-1)-balls[i][2])/(balls[i][9 ]-balls[i][2]))+
				Math.abs(1-(ballZs[i].get(ballZs[i].size()-1)-balls[i][3])/(balls[i][10]-balls[i][3]))
			< 0.01)) {
				balls[i][7] = 1;
			}
			else {
				narrowed = false;
				balls[i][8 ] = ballXs[i].get(ballXs[i].size()-1);
				balls[i][9 ] = ballYs[i].get(ballYs[i].size()-1);
				balls[i][10] = ballZs[i].get(ballZs[i].size()-1);
			}
		}
	}
	// TODO: preserve energy by adjusting velocity here
	// TODO: sort by z-value once colored
	for(int i = 0; i < ball_count; i++) {
		PVector ball = rotate_v((float)balls[i][1], (float)balls[i][2], (float)balls[i][3]);
		circle(ball.x-(float)mx, ball.y-(float)my, (float)balls[i][0]);
		balls[i][1] =  ballXs[i].get( ballXs[i].size()-1);
		balls[i][2] =  ballYs[i].get( ballYs[i].size()-1);
		balls[i][3] =  ballZs[i].get( ballZs[i].size()-1);
		balls[i][4] = ballVXs[i].get(ballVXs[i].size()-1);
		balls[i][5] = ballVYs[i].get(ballVYs[i].size()-1);
		balls[i][6] = ballVZs[i].get(ballVZs[i].size()-1);
	}
	if(three_dimentions) {
		scale((float)(factor/(1-2*win_border)));
		translate(-1+0.15*min(1, width/(float)height), -1+0.15*min(1, height/(float)width));
		scale(0.1);
		strokeWeight(0.1);
		PVector axes[] = new PVector[]{
			new PVector(1,  0, 0),
			new PVector(0, -1, 0),
			new PVector(0,  0, 1)
		};
		if(rotate_z(axes[1]) > rotate_z(axes[0])) { PVector temp = axes[0]; axes[0] = axes[1]; axes[1] = temp; }
		if(rotate_z(axes[2]) > rotate_z(axes[1])) { PVector temp = axes[1]; axes[1] = axes[2]; axes[2] = temp; }
		if(rotate_z(axes[1]) > rotate_z(axes[0])) { PVector temp = axes[0]; axes[0] = axes[1]; axes[1] = temp; }
		for(int i = 0; i < 3; i++) {
			stroke(255*axes[i].x, 255*abs(axes[i].y), 255*axes[i].z);
			PVector axis = rotate_v(axes[i].x, axes[i].y, axes[i].z);
			line(0, 0, axis.x, axis.y);
		}
	}
}
