# Graphics Final Project

I created a 3D animation of a "celestial" watch using a combination of ray tracing, triangle meshes, and 3D rotation matrices. The inspiration for it comes from this watch (https://www.thewatchpages.com/watches/van-cleef-arpels-midnight-planetarium-watch-vcaro4j000/). 

The code for my final project for Computer Graphics (CSCI 371) is in `final_project.cpp`. It needs the `cow.cpp` and `snail.cpp` to run. `cow.cpp` wraps GLFW (a library for OpenGL) and `snail.cpp` is a small linear algebra library. Both were created by Jim Bern who taught CSCI 371. 

Here is are some screenshots of the project: 

<img width="806" alt="graphics_final_project_screenshot" src="https://github.com/cb123450/Graphics/assets/91232059/ea3f8aa4-825e-4b6a-94a3-08e51540b631">

I was able to get shadows to appear on the bezel of the watch, but the animation ran at around 3 fps on my laptop due to the large amount of ray tracing required to do so. Becuase of this, I removed the code for this feature. 
