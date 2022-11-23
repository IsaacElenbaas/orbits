long rand_seed = 1;

// ball measurements are diameter
int   ball_count    = 3;
float ball_min      = 0.1; // relative to smallest window edge
float ball_max_mult = 2.5; // max ball size is this*min
float ball_speed    = 0.1; // max inital ball speed in relative-to-smallest-window-edge/frame
float win_border    = 0.1; // border of screen to ensure is empty, relative to smallest window edge
float max_size      = 40;  // max factor to zoom out by

float gravity = 0.01;

boolean done_init = false;
long cur_rand_seed = rand_seed;
float[][] balls = new float[ball_count][9]; // size, x, y, vx, vy, done for this frame, last substep x, last substep y
ArrayList<Float>[]  ballXs = new ArrayList[ball_count];
ArrayList<Float>[]  ballYs = new ArrayList[ball_count];
ArrayList<Float>[] ballVXs = new ArrayList[ball_count];
ArrayList<Float>[] ballVYs = new ArrayList[ball_count];
ArrayList<Float>[] ballAXs = new ArrayList[ball_count];
ArrayList<Float>[] ballAYs = new ArrayList[ball_count];

/*{{{ float rand()*/
// https://github.com/Qqwy/SimpleRNG
long max_32 = 1L << 32;
float rand() {
	cur_rand_seed ^= (cur_rand_seed << 13)%max_32;
	cur_rand_seed ^= (cur_rand_seed >> 17)%max_32;
	cur_rand_seed ^= (cur_rand_seed <<  5)%max_32;
	return (float)cur_rand_seed;
}
/*}}}*/

float interp(ArrayList<Float> list, float x) {
	if(x == 0) return list.get(1);
	if(x/list.get(0) >= 1) return list.get(list.size()-1);
	float index = (x/list.get(0))*(list.size()-2)+1;
	return (1-index%1)*list.get((int)index)+(index%1)*list.get((int)index+1);
}

/*{{{ void setup()*/
void setup() {
	if(!done_init) {
		for(int i = 0; i < ball_count; i++) {
			ballXs[i]  = new ArrayList<Float>();
			ballYs[i]  = new ArrayList<Float>();
			ballVXs[i] = new ArrayList<Float>();
			ballVYs[i] = new ArrayList<Float>();
			ballAXs[i] = new ArrayList<Float>();
			ballAYs[i] = new ArrayList<Float>();
		}
	}
	for(int i = 0; i < ball_count; i++) {
		for(int j = 0; j <= 10; j++) {
			// size
			balls[i][0] = ball_min+((ball_max_mult-1)*ball_min)*(rand()/max_32);
			// position (edge of circle at -1 to 1)
			balls[i][1] = (1-balls[i][0]/2)*(2*(rand()/max_32)-1);
			balls[i][2] = (1-balls[i][0]/2)*(2*(rand()/max_32)-1);
		}
		// velocity
		balls[i][3] = ball_speed*(rand()/max_32);
		balls[i][4] = sin(rand())*balls[i][3];
		balls[i][3] = cos(rand())*balls[i][3];
	}
	size(1080, 1080, P2D);
}
/*}}}*/

void keyPressed() {
	if(key == ' ') {
		rand_seed++; cur_rand_seed = rand_seed;
		setup();
	}
}

void draw() {
	background(255);
	translate(width/2, height/2);
	strokeWeight(0);
	fill(0);

/*{{{ camera center calculation*/
	// using mean instead of actually centering for a more floaty feel
	float mx = 0, my = 0;
	// mean
	/*for(int i = 0; i < ball_count; i++) {
		mx += balls[i][1];
		my += balls[i][2];
	}
	mx /= ball_count; my /= ball_count;//*/
	// least-squared-error mean, graph x=(some array) and y=sum_0^len of (array[i] - x)^2
	float left_error, right_error;
	float min_x = Float.MAX_VALUE, max_x = -Float.MAX_VALUE;
	float min_y = Float.MAX_VALUE, max_y = -Float.MAX_VALUE;
	for(int i = 0; i < ball_count; i++) {
		min_x = min(min_x, balls[i][1]); max_x = max(max_x, balls[i][1]);
		min_y = min(min_y, balls[i][2]); max_y = max(max_y, balls[i][2]);
	}
	left_error  = 0; for(int i = 0; i < ball_count; i++) left_error  += pow(balls[i][1]-min_x, 2);
	right_error = 0; for(int i = 0; i < ball_count; i++) right_error += pow(balls[i][1]-max_x, 2);
	mx = min_x+((right_error-left_error)/ball_count-pow(max_x-min_x, 2))/(-2*(max_x-min_x));
	left_error  = 0; for(int i = 0; i < ball_count; i++) left_error  += pow(balls[i][2]-min_y, 2);
	right_error = 0; for(int i = 0; i < ball_count; i++) right_error += pow(balls[i][2]-max_y, 2);
	my = min_y+((right_error-left_error)/ball_count-pow(max_y-min_y, 2))/(-2*(max_y-min_y));
/*}}}*/

/*{{{ zoom*/
	float factor = 1;
	for(int i = 0; i < ball_count; i++) {
		balls[i][1] -= mx; balls[i][2] -= my;
		factor = max(factor, max(abs(balls[i][1]), abs(balls[i][2]))+balls[i][0]/2);
	}
	if(factor > max_size) {
		rand_seed++; cur_rand_seed = rand_seed;
		setup();
	}
	scale(min(width, height)/2/(factor/(1-2*win_border)));
/*}}}*/

	for(int i = 0; i < ball_count; i++) balls[i][5] = 0;
	boolean narrowed = false;
	for(int mult_step = 1; !narrowed && mult_step > 0; mult_step = mult_step << 1) {
		for(int i = 0; i < ball_count; i++) {
			if(balls[i][5] != 0) continue;
			 ballXs[i].clear();  ballXs[i].add(0f);  ballXs[i].add(balls[i][1]);
			 ballYs[i].clear();  ballYs[i].add(0f);  ballYs[i].add(balls[i][2]);
			ballVXs[i].clear(); ballVXs[i].add(0f); ballVXs[i].add(balls[i][3]);
			ballVYs[i].clear(); ballVYs[i].add(0f); ballVYs[i].add(balls[i][4]);
			ballAXs[i].clear(); ballAXs[i].add(0f);
			ballAYs[i].clear(); ballAYs[i].add(0f);
		}
		for(float i = 0; i < 1; i += 1f/mult_step) {
			for(int j = 0; j < ball_count; j++) {
				if(balls[j][5] != 0) continue;
				float ax = 0, ay = 0;
				for(int k = 0; k < ball_count; k++) {
					if(k == j) continue;
					// max to prevent infinite acceleration when no collision
					float a = -gravity*balls[k][0]/max(ball_min, pow(interp(ballXs[j], i)-interp(ballXs[k], i), 2)+pow(interp(ballYs[j], i)-interp(ballYs[k], i), 2));
					float theta = atan2(interp(ballYs[j], i)-interp(ballYs[k], i), interp(ballXs[j], i)-interp(ballXs[k], i));
					ax += a*cos(theta); ay += a*sin(theta);
				}
				ballXs[j].add(interp(ballXs[j], i)+(interp(ballVXs[j], i)+ax/mult_step/2)/mult_step); ballXs[j].set(0, i+1f/mult_step);
				ballYs[j].add(interp(ballYs[j], i)+(interp(ballVYs[j], i)+ay/mult_step/2)/mult_step); ballYs[j].set(0, i+1f/mult_step);
				ballAXs[j].add(ax); ballAXs[j].set(0, i);
				ballAYs[j].add(ay); ballAYs[j].set(0, i);
			}
			for(int j = 0; j < ball_count; j++) {
				if(balls[j][5] != 0) continue;
				float ax = 0, ay = 0;
				for(int k = 0; k < ball_count; k++) {
					if(k == j) continue;
					// max to prevent infinite acceleration when no collision
					float a = -gravity*balls[k][0]/max(ball_min, pow(interp(ballXs[j], i+1f/mult_step)-interp(ballXs[k], i+1f/mult_step), 2)+pow(interp(ballYs[j], i+1f/mult_step)-interp(ballYs[k], i+1f/mult_step), 2));
					float theta = atan2(interp(ballYs[j], i+1f/mult_step)-interp(ballYs[k], i+1f/mult_step), interp(ballXs[j], i+1f/mult_step)-interp(ballXs[k], i+1f/mult_step));
					ax += a*cos(theta); ay += a*sin(theta);
				}
				ballVXs[j].add(interp(ballVXs[j], i)+(interp(ballAXs[j], i)+ax)/mult_step/2); ballVXs[j].set(0, i+1f/mult_step);
				ballVYs[j].add(interp(ballVYs[j], i)+(interp(ballAYs[j], i)+ay)/mult_step/2); ballVYs[j].set(0, i+1f/mult_step);
			}
		}
		narrowed = true;
		for(int i = 0; i < ball_count; i++) {
			if(balls[i][5] != 0) continue;
			if(mult_step != 1 && (mult_step >= 64 ||
				abs(1-(ballXs[i].get(ballXs[i].size()-1)-balls[i][1])/(balls[i][6]-balls[i][1]))+
				abs(1-(ballYs[i].get(ballYs[i].size()-1)-balls[i][2])/(balls[i][7]-balls[i][2]))
			< 0.01)) {
				balls[i][5] = 1;
			}
			else {
				narrowed = false;
				balls[i][6] = ballXs[i].get(ballXs[i].size()-1);
				balls[i][7] = ballYs[i].get(ballYs[i].size()-1);
			}
		}
	}
	for(int i = 0; i < ball_count; i++) {
		circle(balls[i][1], balls[i][2], balls[i][0]);
		balls[i][1] =  ballXs[i].get( ballXs[i].size()-1);
		balls[i][2] =  ballYs[i].get( ballYs[i].size()-1);
		balls[i][3] = ballVXs[i].get(ballVXs[i].size()-1);
		balls[i][4] = ballVYs[i].get(ballVYs[i].size()-1);
	}
}
