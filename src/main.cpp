#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <cmath>
#include <SDL2/SDL.h>
#include <bits/stdc++.h>
using namespace std;

//Setting necessary constants

//Gravity in m/s^2
const double GRAVITY = 10;
//Total simulation time
const double SIM_TIME = 10;
const double elasticity = 0.95;

//Interval after which to print environment state
//const double POLL_INTERVAL = 0.2;
//Flag to simulate w.r.t real world clock or simulate as fast as possible
//const bool REALTIME = false;
//Flag to clear screen
//const bool CLEAR_SCREEN = true;
//Class to represent all rigid bodies with finite mass

//Example ball class
class Ball
{
    
    
	public:
        double radius;
        double mass;
        pair<double, double> location;
		pair<double, double> velocity;
		pair<double, double> acceleration;


		Ball(double rad, double mas, double x, double y)
		{
			//Gravity applied by default on all objects
			acceleration = make_pair(0, -1 * GRAVITY);
			velocity = make_pair(0, 0);
			location = make_pair(x, y);
            radius = rad;
            mass = mas;
		}
		pair<double, double> getLocation()
		{
			return location;
		}

		//Functions to update postion based on basic physics principles
		void updateVelocity(double time)
		{
			velocity.first += acceleration.first * time;
			velocity.second += acceleration.second * time;
		}

		void updateLocation(double time)
		{
			location.first += 0.5 * acceleration.first * time * time + velocity.first * time;
			location.second += 0.5 * acceleration.second * time * time + velocity.second * time;
		}

		void updateState(double time)
		{
			updateLocation(time);
			updateVelocity(time);
		}

};


//Class for SDL2 Window graphics
class Framework{
private:
    int height;     // Height of the window
    int width;      // Width of the window
    SDL_Renderer *renderer = NULL;      // Pointer for the renderer
    SDL_Window *window = NULL;      // Pointer for the window
    double totalTime;
    int r=0,g=0,b=255;
    int count=0;
public:
    // Contructor which initialize the parameters.
    Framework(int height_, int width_): height(height_), width(width_){
        SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
        SDL_CreateWindowAndRenderer(width, height, 0, &window, &renderer);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
        SDL_RenderClear(renderer);      // Clear the newly created window
        SDL_RenderPresent(renderer);    // Reflects the changes done in the window.
    }

    // Destructor
    ~Framework(){
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    void draw_circle(int center_x, int center_y, int radius_){
        // Setting the color to be RED with 100% opaque (0% trasparent).
        SDL_SetRenderDrawColor(renderer, r, g, b, 255);

        // Drawing circle
        for(int x=center_x-radius_; x<=center_x+radius_; x++){
            for(int y=center_y-radius_; y<=center_y+radius_; y++){
                if((std::pow(center_y-y,2)+std::pow(center_x-x,2)) <= std::pow(radius_,2)){
                    SDL_RenderDrawPoint(renderer, x, y);
                }
            }
        }
        
    }

    void start_sim(){        
        SDL_Event event;    // Event variable

        auto beforeTime = chrono::steady_clock::now();
        auto curTime = beforeTime;
        totalTime = 0;
        //double interval = POLL_INTERVAL;
        double interval2 = 0;
        double elapsedTime;

        vector<Ball *> entities;

        //Creating the ball
        int n=5;
        Ball *ball[n];
        ball[0]=new Ball(0.5,1,0,5);
        ball[1]=new Ball(0.5,2,+2,3.5);
        ball[2]=new Ball(0.5,1.3,-1,4);
        ball[3]=new Ball(0.2,0.7,-1,6);
        ball[4]=new Ball(0.15,0.2,+1,3);

        ball[0]->velocity.first=0;
        ball[1]->velocity.first=-2;
        ball[2]->velocity.first=3;
        ball[3]->velocity.first=-1;
        ball[4]->velocity.first=1;
        for(int i=0;i<n;i++){
            entities.push_back(ball[i]);
        }
        //int count=0;
        while(!(event.type == SDL_QUIT)){
            //Set background to black
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
            SDL_RenderClear(renderer);
            // SDL_Delay(1);
            curTime = chrono::steady_clock::now();
			elapsedTime = chrono::duration<double, micro>(curTime - beforeTime).count() / 1e6;
			beforeTime = curTime;
            totalTime += elapsedTime;
            interval2+=elapsedTime;
            int split=100;
            //Updating the state of all entities
            for(int j=0;j<split;j++){
                for (int i = 0; i < entities.size(); i++)
                {
                    entities[i]->updateState(elapsedTime/split);
                }

                //Collision detection
                checkCollision(entities);
            }



            for (int i = 0; i < entities.size(); i++)
            {
                int y=(int)(entities[i]->getLocation().second*100);
                int x=(int)(entities[i]->getLocation().first*100);

                //Draw circle with respective coordinates
                draw_circle(width/2+x, height-y, entities[i]->radius*100);
            }            
            //Render screen
            SDL_RenderPresent(renderer);

            //End simulation when time complete
            if (totalTime >= SIM_TIME){
                cout << "FPS: " << count/SIM_TIME << endl;
                break;
            }
            count++;
            // if(interval2>interval){
            //     r=random()%256;
            //     g=random()%256;
            //     b=random()%256;
            //     interval2-=interval;
            // }
            //Check for window close action    
            SDL_PollEvent(&event);
           
        }
        
    }

    //Function to detect and handle collisions in the environment
    //Currently there is a ground plane at y postion 0, and the coefficient of elasticity is assumed to be unity, so the ball bounces back with the same speed as it hits the ground
    void checkCollision(vector<Ball *> entities)
    {
        for (int i = 0; i < entities.size(); i++)
        {
            //Get radius of first ball
            double r = entities[i] -> radius;

            //Iterate through all pairs
            for(int j = 0; j< entities.size(); j++){
                if(i==j) continue;
                //Masses
                double m1=entities[i]->mass, m2=entities[j]->mass;
                //Radius of second ball
                double r2 = entities[j] -> radius;

                //Distance btw centers
                double dis = pow(pow(entities[j]->location.second-entities[i]->location.second,2)+pow(entities[j]->location.first-entities[i]->location.first,2),0.5);
                
                //If colliding
                if(dis<r+r2){
                    //Newtons equations
                    double x1 = entities[i]->location.first;
                    double x2 = entities[j]->location.first;
                    double y1 = entities[i]->location.second;
                    double y2 = entities[j]->location.second;
                    double xv1 = entities[i]->velocity.first;
                    double xv2 = entities[j]->velocity.first;
                    double yv1 = entities[i]->velocity.second;
                    double yv2 = entities[j]->velocity.second;

                    // Calculate displacement required
                    double fOverlap = 0.5 * (dis - r - r2);

                    // Displace Current Ball away from collision
                    entities[i]->location.first -= fOverlap * (x1 - x2) / dis;
                    entities[i]->location.second -= fOverlap * (y1 - y2) / dis;

                    // Displace Target Ball away from collision
                    entities[j]->location.first += fOverlap * (x1 - x2) / dis;
                    entities[j]->location.second += fOverlap * (y1 - y2) / dis;

                    dis = pow(pow(entities[j]->location.second-entities[i]->location.second,2)+pow(entities[j]->location.first-entities[i]->location.first,2),0.5);

                    //Normal
                    double nx = (entities[j]->location.first-entities[i]->location.first)/dis;
                    double ny = (entities[j]->location.second-entities[i]->location.second)/dis;
                    //Tangent
                    double tx = -ny;
                    double ty = nx;

                    //Dot product tangent
                    double dpTan1 = xv1*tx + yv1*ty;
                    double dpTan2 = xv2*tx + yv2*ty;

                    //Dot product normal
                    double dpNorm1 = xv1*nx + yv1*ny;
                    double dpNorm2 = xv2*nx + yv2*ny;

                    double mo1 = (dpNorm1 * (m1 - elasticity*m2) + (1+elasticity)* m2 * dpNorm2) / (m1+m2);
			        double mo2 = (dpNorm2 * (m2 - elasticity*m1) + (1+elasticity) * m1 * dpNorm1) / (m1+m2);


                    entities[i]->velocity.first = tx*dpTan1 + mo1*nx;
                    entities[i]->velocity.second = ty*dpTan1 + mo1*ny;
                    entities[j]->velocity.first = tx*dpTan2 + mo2*nx;
                    entities[j]->velocity.second = ty*dpTan2 + mo2*ny;
                    
                }
            }
            //Floor collision check
            if (entities[i]->location.second < r)
            {
                double a = entities[i]->acceleration.second;
                double v = entities[i]->velocity.second;
                double u = -pow(v*v-2*a*(entities[i]->location.second - r),0.5);
                double t = (v-u)/a;
                entities[i]->location.second=max(r,r+(-1*elasticity*u)*t+(0.5*a*t*t));
                entities[i]->velocity.second=max(0.0,(-1*elasticity*u)+a*t);
            }

            //Upper wall
            if (entities[i]->location.second > height/100.0-r && entities[i]->velocity.second >= 0)
            {
                entities[i]->velocity.second *= -1 * elasticity;
            }

            //Left wall
            if (entities[i]->location.first <= -1 * (width/200.0 - r) && entities[i]->velocity.first <= 0)
            {
                entities[i]->velocity.first *= -1 * elasticity;
            }

            //Right wall
            if (entities[i]->location.first >= (width/200.0 - r) && entities[i]->velocity.first >= 0)
            {
                entities[i]->velocity.first *= -1 * elasticity;
            }
        }
    }

};


signed main()
{
	ios::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);

    // Creating the object by passing Height and Width value.
    Framework fw(750, 750);

    // Starting the animation
    fw.start_sim();

	return 0;
}