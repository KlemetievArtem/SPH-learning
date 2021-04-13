#include "Application.h"

int main(void){

	   


	Application application("tutorial", 1920, 1080, 4, 6, false);

	//MAIN LOOP
	while (!application.getWindowShouldClouse()) {
		//UPDATE INPUT
		application.update();
		//application.SetEventCallback();


		//œ¿–¿ÀÀ≈À»“‹
		application.render();
	}


	return 0;
}