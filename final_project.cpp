#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"
#include <string> 

char *final_frag = R""""(
    
    //SAVESAVE
    #version 330 core

    // https://iquilezles.org/articles/distfunctions/
    float dot2(vec2 v) { return dot(v,v); }
    float dot2(vec3 v) { return dot(v,v); }
    float ndot(vec2 a, vec2 b) { return a.x*b.x - a.y*b.y; }

    // for computing ray directions
    uniform vec3 x_renderer;
    uniform vec3 y_renderer;
    uniform vec3 z_renderer;
    uniform vec3 o_renderer;
    uniform float renderer_angle_of_view;
    uniform vec2 iResolution;
    
    uniform vec3 planet_positions[9];
    uniform vec3 planet_colors[9];
    uniform vec3 planet_radii[3];
    uniform int speed;

    uniform vec3 planet_change_theta[4];

    uniform vec3 asteroid_positions[48];
    uniform float asteroid_radius;
    uniform bool blinnphong; 
    uniform bool shadow_on_ring; 
    uniform int intensity; 

    uniform float time; // for time-varying distance fields

    out vec4 fragColor;
    
    // begin https://iquilezles.org/articles/distfunctions/
    float sdSphere(vec3 p, float r) {
        return length(p) - r;
    }

    float sdTorus(vec3 p, vec2 t)
    {
        vec2 q = vec2(length(p.xz)-t.x,p.y);
        return length(q)-t.y;
    }

    struct ColorGlow{
        vec4 col;
        float glow;
    };

    vec3 light_color = vec3(1., 1., 1.); //white light
    vec3 light_position = planet_positions[0]; //the sun is the light
    vec3 sun_position = light_position; 
    float sun_radius = planet_radii[0][0];

    //For shading!
    vec3 normal_to_sphere(vec3 p, float sphere_radius){
        //see https://michaelwalczyk.com/blog-ray-marching.html
        const vec3 epsilon = vec3(0.001, 0.0, 0.0);

        //get gradient at point 
        //swizzling
        float x_gradient = sdSphere(p+epsilon.xyy, sphere_radius) - sdSphere(p-epsilon.xyy, sphere_radius);
        float y_gradient = sdSphere(p+epsilon.yxy, sphere_radius) - sdSphere(p-epsilon.yxy, sphere_radius);
        float z_gradient = sdSphere(p+epsilon.yyx, sphere_radius) - sdSphere(p-epsilon.yyx, sphere_radius);

        vec3 normal = vec3(x_gradient, y_gradient, z_gradient);

        return normalize(normal);
    }

    //For shading!
    vec3 normal_to_torus(vec3 p, vec2 t){
        /*see https://www.shadertoy.com/view/ll33Wn but it is essentially the same as finding 
        the normal to a sphere; just get the gradient to a point on the sphere */
        const vec3 epsilon = vec3(0.001, 0.0, 0.0);

        //get gradient at point 
        //swizzling
        float x_gradient = sdTorus(p+epsilon.xyy, t) - sdTorus(p-epsilon.xyy, t);
        float y_gradient = sdTorus(p+epsilon.yxy, t) - sdTorus(p-epsilon.yxy, t);
        float z_gradient = sdTorus(p+epsilon.yyx, t) - sdTorus(p-epsilon.yyx, t);

        vec3 normal = vec3(x_gradient, y_gradient, z_gradient);

        return normalize(normal);

    }
    
    //for shadows, march from the planet to the sun
    //1 if hit a planet or an asteroid and 0 if do not hit a planet
    //SHOULD BE USING DISTANCE-AIDED RAY MARCHING AND INCREMENT
    //'total_distance_marched' by the distance to the closest object each time

    int march_from_planet_to_sun(vec3 start, vec3 dir, int index){

        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = .001;
        float MAX_MARCH_DISTANCE = length(dir);

        // next   -- current position along ray
        // total_distance_marched   -- distance marched along ray
        // dir -- distance in direction of planet from the sun          
        // f   -- distance to surface       
        
        float total_distance_marched = 0.0;
        int step = 0;

        while (step++ < MAX_STEPS && total_distance_marched < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 next = start + total_distance_marched * dir;

            // compute distance to implicit surface
            float f = MAX_MARCH_DISTANCE;


            //ignore the sun (current planet is ignored because of offset)
            //consider planets
            for (int planet = index; planet > 0; planet--){
                int i = planet/3;
                int k = planet % 3; 

                mat3 planet_rotation_matrix = mat3(cos(planet_change_theta[i][k]*time), -sin(planet_change_theta[i][k]*time), 0,
                                        sin(planet_change_theta[i][k]*time), cos(planet_change_theta[i][k]*time), 0,
                                        0, 0, 1);

                vec3 planet_pos = planet_rotation_matrix*planet_positions[planet];

                float sphere_radius = planet_radii[i][k];

                float distance_to_sphere = sdSphere(next - planet_pos, sphere_radius);
                f = min(f, distance_to_sphere);

                if (f < HIT_TOLERANCE && length(start- 0) > length(planet_pos-0)) { // hit!
                    //ensures that the planet casting the shadow is closer to the sun than the planet getting the shadow cast on 
                    return 1;
                }
            }

            //consider asteroids
            for (int i = 0; i < 48; i++){
                mat3 asteroid_rotation_matrix = mat3(cos(planet_change_theta[3][0]*time), -sin(planet_change_theta[3][0]*time), 0,
                        sin(planet_change_theta[3][0]*time), cos(planet_change_theta[3][0]*time), 0,
                        0, 0, 1);

                vec3 asteroid_pos = asteroid_rotation_matrix*asteroid_positions[i];

                float sphere_radius = asteroid_radius;

                float distance_to_sphere = sdSphere(next - asteroid_pos, sphere_radius);
                f = min(f, distance_to_sphere);

                if (f < HIT_TOLERANCE && length(start- 0) > length(asteroid_pos-0)) { // hit!
                    //ensures that the planet casting the shadow is closer to the sun than the planet getting the shadow cast on 
                    return 1;
                }
            }
                
            // NOTE if you're getting weird "overstepping" artifacts
            // (weird missing slices in the geometry)               
            // a (hacky) solution is to replace t += f; with e.g.   
            //total_distance_marched += min(f, .01);
            total_distance_marched += f / 15.;
        }
        //did not hit
        return 0;
    }

    
    //returns 0 if asteroid is in the shadow of a planet and returns 1 otherwise
    int asteroid_march_to_sun(vec3 start, vec3 dir, vec3 asteroid_pos){

        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = .001;
        float MAX_MARCH_DISTANCE = length(dir);

        // next   -- current position along ray
        // total_distance_marched   -- distance marched along ray
        // dir -- distance in direction of planet from the sun          
        // f   -- distance to surface       
        
        float total_distance_marched = 0.0;
        int step = 0;

        while (step++ < MAX_STEPS && total_distance_marched < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 next = start + total_distance_marched * dir;

            // compute distance to implicit surface
            float f = MAX_MARCH_DISTANCE;

            //ignore the sun (current asteroid is ignored because of offset)
            //consider PLANETS WIHTIN RANGE MERCURY TO MARS
            for (int planet = 4; planet > 0; planet--){
                int i = planet/3;
                int k = planet % 3; 

                mat3 planet_rotation_matrix = mat3(cos(planet_change_theta[i][k]*time), -sin(planet_change_theta[i][k]*time), 0,
                                        sin(planet_change_theta[i][k]*time), cos(planet_change_theta[i][k]*time), 0,
                                        0, 0, 1);

                vec3 planet_pos = planet_rotation_matrix*planet_positions[planet];

                float sphere_radius = planet_radii[i][k];

                float distance_to_sphere = sdSphere(next - planet_pos, sphere_radius);
                f = min(f, distance_to_sphere);

                if (f < HIT_TOLERANCE && length(start- 0) > length(planet_pos-0)) { // hit!
                    //ensures that the planet casting the shadow is closer to the sun than the planet getting the shadow cast on 
                    return 1;
                }
            }

            /*
            //commented out because asteroids will never block each others path to the sun based on how I placed them; this speeds it up by 5 fps about
            //consider asteroids
            for (int i = 0; i < 48; i++){
                mat3 asteroid_rotation_matrix = mat3(cos(time*planet_change_theta[3][0]), -sin(time*planet_change_theta[3][0]), 0,
                        sin(time*planet_change_theta[3][0]), cos(time*planet_change_theta[3][0]), 0,
                        0, 0, 1);

                vec3 asteroid_pos = asteroid_rotation_matrix*asteroid_positions[i];

                float sphere_radius = asteroid_radius;

                float distance_to_sphere = sdSphere(next - asteroid_pos, sphere_radius);
                f = min(f, distance_to_sphere);

                if (f < .001*HIT_TOLERANCE && length(start- 0) > length(asteroid_pos-0)) { // hit!
                    //ensures that the planet casting the shadow is closer to the sun than the planet getting the shadow cast on 
                    return 1;
                }
            }
            */
                
            // NOTE if you're getting weird "overstepping" artifacts
            // (weird missing slices in the geometry)               
            // a (hacky) solution is to replace t += f; with e.g.   
            // t += min(f, .5);
            total_distance_marched += f / 15.;
        }
        //did not hit
        return 0;
    }
    
    

    vec4 march(vec3 o, vec3 dir) {
        // I used this for reference: https://michaelwalczyk.com/blog-ray-marching.html

        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = 0.001;
        const float MAX_MARCH_DISTANCE = 100.0;

        // p   -- current position along ray
        // o   -- camera origin             
        // total_distance_marched   -- distance marched along ray
        // dir -- camera direction          
        // f   -- distance to surface       

        float total_distance_marched = 0.0;
        int step = 0;
        
        //the sun is the source of light

        while (step++ < MAX_STEPS && total_distance_marched < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 p = o + total_distance_marched * dir;

            // compute distance to implicit surface
            float f = MAX_MARCH_DISTANCE; {
                
                //PLANETS
                {   
                    
                    //make the sun glow
                    vec3 sun_position = planet_positions[0];
                    float distance_to_sun = sdSphere(p - sun_position, sun_radius);
                    f = min(f, distance_to_sun);

                    if (f < HIT_TOLERANCE){
                        return vec4(planet_colors[0], 1.0);
                    }

                    for (int planet = 1; planet < 9; planet++){
                        int i = planet/3;
                        int k = planet % 3;
                        
                        vec3 planet_pos = planet_positions[planet];
                        if (planet > 0){
                            mat3 planet_rotation_matrix = mat3(cos(planet_change_theta[i][k]*time), -sin(planet_change_theta[i][k]*time), 0,
                                    sin(planet_change_theta[i][k]*time), cos(planet_change_theta[i][k]*time), 0,
                                    0, 0, 1);

                            planet_pos = planet_rotation_matrix*planet_positions[planet];
                        }

                        vec3 sphere_position = planet_positions[planet];
                        float sphere_radius = planet_radii[i][k];
                        float distance_to_sphere = sdSphere(p - planet_pos, sphere_radius);
                        f = min(f, distance_to_sphere);

                        if (f < HIT_TOLERANCE) { // hit!
                            //ambient lighting
                            float ambientStrength = .25;
                            vec3 ambient = ambientStrength*light_color;

                            //get normal to point on sphere; normals point outwards
                            vec3 normal = normal_to_sphere(p-planet_pos, planet_radii[i][k]);

                            //march from the planet to the sun
                            
                            vec3 start = (planet_positions[0]-p)/1000.+p;
                            if (march_from_planet_to_sun(start, planet_positions[0]-p, planet) == 0){
                                //ONLY ADD DIFFUSE AND SPECULAR LIGHTING TO A PLALENT IF IT IS NOT IN THE SHADOW OF ANOTHER PLANET OR ASTEROID

                                //the sun is the source of light
                                vec3 light_dir = normalize(light_position-p);

                                // Phong lighting: see https://learnopengl.com/Lighting/Basic-Lighting 

                                //diffuse lighting
                                float diff = max(dot(normal,light_dir), 0.0);
                                vec3 diffuse = diff*light_color;

                                //specular lighting
                                float specular_strength = .7;

                                vec3 view_dir = normalize(o_renderer - p);
                                vec3 reflect_dir = reflect(-light_dir, normal);
                                int shininess = 32;
                                vec3 halfway_dir;

                                float spec = pow(max(dot(view_dir, reflect_dir), 0.0), shininess);
                                vec3 specular = specular_strength*spec*light_color;

                                
                                if (blinnphong){
                                    view_dir = normalize(o_renderer - p);

                                    vec3 halfway_dir = normalize(light_dir+view_dir);
                                    spec = pow(max(dot(normal, halfway_dir), 0.0), shininess);
                                    vec3 specular = specular_strength*spec*light_color;
                                }

                                if (planet > 0){
            
                                    vec3 result = (ambient+diffuse+specular)*planet_colors[planet];
                                    return vec4(result, 1.);
                                }
                                //the sun glows
                                
                            }
                            else{
                                ambient += vec3(-.1, -.1, 0.); 
                                //this makes the background of planets less yellow to counterbalance the yellow glow from the sun, 
                                //which hits every pixel but should not light the back of any planets/asteroids

                                vec3 result = (ambient)*planet_colors[planet];
                                return vec4(result, 1.);
                            }
                        }
                    }
                }

                //ASTEROIDS
                {   
                    for (int i = 0; i < 48; i++){
                        mat3 asteroid_rotation_matrix = mat3(cos(time*planet_change_theta[3][0]), -sin(time*planet_change_theta[3][0]), 0,
                                sin(time*planet_change_theta[3][0]), cos(time*planet_change_theta[3][0]), 0,
                                0, 0, 1);

                        vec3 asteroid_pos = asteroid_rotation_matrix*asteroid_positions[i];

                        //gray
                        vec3 asteroid_color = vec3(.5, .5, .5);

                        //float asteroid_radius = .1;
                        float distance_to_sphere = sdSphere(p - asteroid_pos, asteroid_radius);
                        f = min(f, distance_to_sphere);

                        if (f < HIT_TOLERANCE) { // hit!
                            //get normal to point on sphere
                            vec3 normal = normal_to_sphere(p-asteroid_pos, asteroid_radius);

                            //the sun is the source of light
                            vec3 light_dir = normalize(light_position-p);

                            vec3 result;

                            // Phong lighting: see https://learnopengl.com/Lighting/Basic-Lighting 
                            //ambient lighting
                            float ambientStrength = .25;
                            vec3 ambient = ambientStrength*light_color;

                            //diffuse lighting
                            float diff = max(dot(normal,light_dir), 0.0);
                            vec3 diffuse = diff*light_color;
                            
                            //specular lighting
                            float specular_strength = .7;

                            vec3 view_dir = normalize(o_renderer - p);
                            vec3 reflect_dir = reflect(-light_dir, normal);
                            int shininess = 32;
                            vec3 halfway_dir;

                            float spec = pow(max(dot(view_dir, reflect_dir), 0.0), shininess);
                            vec3 specular = specular_strength*spec*light_color;
                            
                            if (blinnphong){
                                view_dir = normalize(o_renderer-p);

                                vec3 halfway_dir = normalize(light_dir+view_dir);
                                spec = pow(max(dot(normal, halfway_dir), 0.0), shininess);
                                vec3 specular = specular_strength*spec*light_color;
                            }


                            vec3 start = (planet_positions[0]-p)/1000.+p;

                            //check if an asteroid is in the shadow of another asteroid or a planet with 1 <= index <=4  (MARS, EARTH, VENUS, MERCURY)
                            if (asteroid_march_to_sun(start, planet_positions[0]-start, asteroid_pos) == 0){
                                result = (ambient+diffuse+specular)*asteroid_color;
                            }
                            else{
                                ambient += vec3(-.1, -.1, 0.); 
                                result = (ambient)*asteroid_color;
                            }
                            return vec4(result, 1.);
                        }
                    }
                }

                //WATCH RING
                {
                    //rotate torus to align it with rest of watch
                    mat3 rot_x = mat3(1, 0, 0,
                                        0, cos(3.14/2.0), -sin(3.14/2.0),
                                        0, sin(3.14/2.0), cos(3.14/2.0) );

                    //for metal part of watch face
                    vec3 new_torus_position = rot_x*vec3(0.0, 0.0, -.5);
                    vec3 new_p = rot_x*p;

                    float torus_major_radius = 19.0;
                    float torus_minor_radius = 1.0;
                    vec3 torus_color = vec3(.5, .5, .5);
                    vec2 torus_radii = vec2(torus_major_radius, torus_minor_radius);
                    float distance_to_torus = sdTorus(new_p - new_torus_position, torus_radii);
                    f = min(f, distance_to_torus);
                    
                    
                    if (f < HIT_TOLERANCE){
                        return vec4(torus_color, 1.);
                    }

                    /*
                    //RUNS ON 5 FPS IF SHADOWS ARE ON RING
                    if (!shadow_on_ring){
                        if (f < HIT_TOLERANCE){
                            return vec4(torus_color, 1.);
                        }
                    }
                    if (f < HIT_TOLERANCE){
                        if (!shadow_on_ring){
                            return vec4(torus_color, 1.);
                        }
                        else{
                            //ambient lighting
                            float ambientStrength = .25;
                            vec3 ambient = ambientStrength*light_color;

                            //get normal to point on sphere; normals point outwards
                            vec3 normal = normal_to_torus(new_p-new_torus_position, torus_radii);

                            vec3 start = (planet_positions[0]-p)/1000.+p;
                            if (march_from_planet_to_sun(start, planet_positions[0]-p, 9) == 0){
                                //ONLY ADD DIFFUSE AND SPECULAR LIGHTING TO THE WATCH OUTLINE IF IT IS NOT IN THE SHADOW OF ANOTHER PLANET OR ASTEROID

                                //the sun is the source of light
                                vec3 light_dir = normalize(light_position-p);

                                // Phong lighting: see https://learnopengl.com/Lighting/Basic-Lighting 

                                //diffuse lighting
                                float diff = max(dot(normal, light_dir), 0.0);
                                vec3 diffuse = diff*light_color;

                                //specular lighting
                                float specular_strength = .7;

                                vec3 view_dir = normalize(o_renderer - p);
                                vec3 reflect_dir = reflect(-light_dir, normal);
                                int shininess = 32;
                                vec3 halfway_dir;

                                float spec = pow(max(dot(view_dir, reflect_dir), 0.0), shininess);
                                vec3 specular = specular_strength*spec*light_color;
                                
                                if (blinnphong){
                                    view_dir = normalize(o_renderer - p);

                                    vec3 halfway_dir = normalize(light_dir+view_dir);
                                    spec = pow(max(dot(normal, halfway_dir), 0.0), shininess);
                                    vec3 specular = specular_strength*spec*light_color;
                                }

                                vec3 result = (ambient+diffuse+specular)*torus_color;
                                return vec4(result, 1.);
                                
                            }
                            else{
                                //this makes the background of planets less yellow to counterbalance the yellow glow from the sun, 
                                //which hits every pixel but should not light the back of any planets/asteroids

                                vec3 result = (ambient)*torus_color;
                                return vec4(result, 1.);
                            }
                        }
                    }
                    */
                    
                }
            }

            // NOTE if you're getting weird "overstepping" artifacts
            // (weird missing slices in the geometry)               
            // a (hacky) solution is to replace t += f; with e.g.   
            // total_distance_marched += min(f, .5);
            total_distance_marched += f / 2.;
        }
        return vec4(0.0);
    }

    float accumulateGlow(vec3 ray_origin, vec3 dir){
        //I based my 3D glow implementation on this example on shadertoy https://www.shadertoy.com/view/7stGWj 
        float glow = 0;
        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = 0.001;
        const float MAX_MARCH_DISTANCE = 100.0;

        float total_distance_marched = 0.0;
        int step = 0;

        while (step++ < MAX_STEPS && total_distance_marched < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 current_position_on_ray = ray_origin + total_distance_marched * dir;

            // compute distance to implicit surface
            float closest_distance_to_sun = MAX_MARCH_DISTANCE; {

                //accumulate glow; rays closer to the sun get more glow 
                float distance_to_sun = sdSphere(current_position_on_ray - sun_position, sun_radius);
                closest_distance_to_sun = min(closest_distance_to_sun, distance_to_sun);
                glow += pow(1e-3 / max(closest_distance_to_sun, HIT_TOLERANCE), intensity/100.);
                
                if (closest_distance_to_sun < HIT_TOLERANCE){
                    return -1;
                }
            }
            total_distance_marched += closest_distance_to_sun /2.;
        }
        return glow;
    }
    
    

    void main() {
        vec3 o = o_renderer;
        vec3 dir; {
            // NOTE assume unit distance to film plane
            vec2 ds; { // [-R, R]
                float theta = renderer_angle_of_view / 2;
                float _R = tan(theta);
                ds = gl_FragCoord.xy;
                ds -= vec2(iResolution.x / 2, iResolution.y / 2);
                ds *= _R * 2. / iResolution.y;
            }
            dir = -z_renderer + ds.x * x_renderer + ds.y * y_renderer;
        }
        
        vec4 col = march(o, dir);

        //call separate glow function to make code easier to follow and 
        float glow = accumulateGlow(o, dir);

        //the clamp prevents the white line at the radius of the sun (the part with most glow) from being white
        glow = clamp(glow, -1., .99999); 

        if (glow == -1){
            fragColor = col;
        }
        else{
            vec4 glow_color = vec4(planet_colors[0], 1.0); //yellow
            //vec4 glow_component = glow*glow_color;
            //glow is never accumulated
            fragColor = col;

            fragColor = col + glow*glow_color;
        }
    }
)"""";
    

void finalproject() {
    init();

    //BEGIN SHADER SETUP
    // mesh
    int num_vertices = 4;
    int num_triangles = 2;
    vec3 vertex_positions[] = { { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 }, { -1, 1, 0 } };
    int3 triangle_indices[] = { { 0, 1, 2 }, { 0, 2, 3 } };

    // shaders
    char *vert = R""""(
        #version 330 core
        layout (location = 0) in vec3 _p_model;
        void main() {
            gl_Position = vec4(_p_model, 1);
        }
    )"""";

    int shader_program = shader_build_program(vert, final_frag);
    ASSERT(shader_program);

    // misc opengl
    GLuint VAO, VBO, EBO; {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
    }
    //END SHADER SETUP    

    //BEGIN CAMERA SETUP
    // misc cow
    bool playing = false;
    bool shadow_on_ring = false;
    double time = 0;

    double screen_height_world = 18.0;
    double angle_of_view = 3.14/4;

    double sphere_radius = (screen_height_world / 2) / tan(angle_of_view / 2); 

    double theta = 0.0; // yaw
    double phi = 0.0; // pitch

    Camera3D renderer = {sphere_radius, angle_of_view, theta, phi, 0.0, 0.0};
    //END CAMERA SETUP


    //BEGIN PLANET SETUP
    int const num_planets = 9;
    //initialize planets
    //sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune
    vec3* planet_positions = (vec3*) calloc(num_planets, sizeof(vec3));
    vec3* planet_radii = (vec3*) calloc(num_planets/3, sizeof(vec3));
    vec3* planet_colors = (vec3*) calloc(num_planets, sizeof(vec3)); 
    
    planet_colors[0] = monokai.yellow;
    planet_colors[1] = monokai.red;
    planet_colors[2] = V3(.5, .5, .3);
    planet_colors[3] = monokai.blue;
    planet_colors[4] = monokai.red;
    planet_colors[5] = monokai.orange;
    planet_colors[6] = monokai.green;
    planet_colors[7] = V3(.2, .3, .54);
    planet_colors[8] = V3(.0, .3, .54);
    
    planet_radii[0] = V3(1.5, .20, .35)/2.0;
    planet_radii[1] = V3(.40, .35, 2.0)/2.0;
    planet_radii[2] = V3(1.50, .40, .45)/2.0;

    float distance_from_origin = 0.0;
    int index = 0;
    double inner_aster_radius = 0;
    double outer_aster_radius = 0;

    for (int i = 0; i < num_planets/3; i++){
        for (int k = 0; k < num_planets/3; k++){
            planet_positions[index] = V3(distance_from_origin, 0.0, 0.0);
            distance_from_origin += planet_radii[i][k] + 1.25;
            if (index == 3){
                inner_aster_radius = distance_from_origin;
                distance_from_origin += 2.5;
                outer_aster_radius = distance_from_origin;
                distance_from_origin += .5;
            }
            index += 1;
        }
    }


    //do not rotate the sun
    vec3* planet_change_theta_orig = (vec3*) calloc(1+ num_planets/3, sizeof(vec3));

    planet_change_theta_orig[0] = V3(0.0, .080, .045);
    planet_change_theta_orig[1] = V3(.025, .090, .050);
    planet_change_theta_orig[2] = V3(.030, .040, .035);
    planet_change_theta_orig[3] = V3(.030, 0., 0.);

    /*
    0: SUN
    1: MERCURY
    2: VENUS
    3: EARTH
    4: MARS
    5: JUPITER
    6: SATURN
    7: URANUS
    8: NEPTUNE
    */

    //END PLANET SETUP

    //ASTEROIDS BETWEEN PLANET AT INDEX 4 (MARS) AND PLANET AT INDEX 5 (JUPITER)
    //START ASTEROID BELT
    int num_asteroids = 3*16;
    vec4* asteroid_positions_hom = (vec4*) calloc(num_asteroids, sizeof(vec4));


    float asteroid_radius = .1;

    //3 concentric circles of asteroids
    double separation = (outer_aster_radius-inner_aster_radius)/3.0;

    double change_theta = 2*3.1415926/16;
    double start_theta[3] = {0, 4/3.14, 2/3.14}; 
    for (int i = 0; i < 16; i++){
        double start = inner_aster_radius; 
        for (int k = 0; k < 3; k++){
            asteroid_positions_hom[3*i+k] = RotationZ(start_theta[k])*V4(start, 0.0, 0.0, 1.0);
            //need homogenous coordinates to rotate
            start += separation;

            //rotate a little so 3 asteroids are not all on the same line
            start_theta[k] += change_theta; 
        }
    }

    vec3* asteroid_positions = (vec3*) calloc(num_asteroids, sizeof(vec3));
    for (int i = 0; i < num_asteroids; i++){
        vec4 hom_asteroid_i = asteroid_positions_hom[i];
        asteroid_positions[i] = V3(hom_asteroid_i.x, hom_asteroid_i.y, hom_asteroid_i.z);
    }

    //END ASTEROID BELT


    //BEGIN CONSTELLATIONS
    int big_dipper_vert = 7;
    vec4* big_dipper = (vec4*) calloc(big_dipper_vert*2, sizeof(vec4));
    big_dipper[0] = V4(0.0, 1.0, 0, 1.0); //point that is at the intersection between the handle and the spoon
    big_dipper[1] = V4(2.0, 1.0, 0, 1.0);
    big_dipper[2] = V4(2.0, 1.0, 0, 1.0);
    big_dipper[3] = V4(1.8, 0.0, 0, 1.0);
    big_dipper[4] = V4(1.8, 0.0, 0, 1.0);
    big_dipper[5] = V4(0.2, 0.0, 0, 1.0);
    big_dipper[6] = V4(0.2, 0.0, 0, 1.0);
    big_dipper[7] = V4(0.0, 1.0, 0, 1.0);
    //end spoon
    //begin handle 
    big_dipper[8] = V4(0.0, 1.0, 0, 1.0);
    big_dipper[9] = V4(-.5, 1.5, 0, 1.0);
    big_dipper[10] = V4(-.5, 1.5, 0, 1.0);
    big_dipper[11] = V4(-2.0, 2.5, 0, 1.0);
    big_dipper[12] = V4(-2.0, 2.5, 0, 1.0);
    big_dipper[13] = V4(-4.0, 2.0, 0, 1.0);

    for (int i = 0; i < big_dipper_vert*2; i++){
        big_dipper[i] = Translation(3, -3, -2.5)*Scaling(.3, .3, .3)*big_dipper[i];
    }

    
    //a triangle soup mesh and organize counterclockwise
    const int big_dipper_mesh_vertices = 6;
    vec4 big_dipper_triangle_mesh[big_dipper_mesh_vertices] = {
                                big_dipper[0], big_dipper[3], big_dipper[1],
                                big_dipper[0], big_dipper[5], big_dipper[3]
                                };
    

    int ursa_major_vert = 16; 
    vec4* ursa_major = (vec4*) calloc(ursa_major_vert*2, sizeof(vec4));

    ursa_major[0] = V4(0.0, 0.0, 0.0, 1.0);
    ursa_major[1] = V4(0.9, -.11, 0.0, 1.0);

    ursa_major[2] = V4(0.9, -.11, 0.0, 1.0);
    ursa_major[3] = V4(2.1, -1.7, 0.0, 1.0);

    ursa_major[4] = V4(2.1, -1.7, 0.0, 1.0);
    ursa_major[5] = V4(5.7, -2.0, 0.0, 1.0);

    ursa_major[6] = V4(5.7, -2.0, 0.0, 1.0);
    ursa_major[7] = V4(9.6, -2.3, 0.0, 1.0);

    ursa_major[8] = V4(9.6, -2.3, 0.0, 1.0);
    ursa_major[9] = V4(12.0, -3.0, 0.0, 1.0);

    ursa_major[10] = V4(12.0, -3.0, 0.0, 1.0);
    ursa_major[11] = V4(9.2, -4.1, 0.0, 1.0);

    ursa_major[12] = V4(9.2, -4.1, 0.0, 1.0);
    ursa_major[13] = V4(5.2, -4.2, 0.0, 1.0);

    ursa_major[14] = V4(5.2, -4.2, 0.0, 1.0);
    ursa_major[15] = V4(2.6, -4.0, 0.0, 1.0);

    ursa_major[16] = V4(2.6, -4.0, 0.0, 1.0);
    ursa_major[17] = V4(0.0, -6.9, 0.0, 1.0);

    ursa_major[18] = V4(2.6, -4.0, 0.0, 1.0);
    ursa_major[19] = V4(3.2, -5.4, 0.0, 1.0);

    ursa_major[20] = V4(3.2, -5.4, 0.0, 1.0);
    ursa_major[21] = V4(3.9, -8.0, 0.0, 1.0);

    ursa_major[22] = V4(3.2, -5.4, 0.0, 1.0);
    ursa_major[23] = V4(4.4, -7.5, 0.0, 1.0);

    ursa_major[24] = V4(9.2, -4.1, 0.0, 1.0);
    ursa_major[25] = V4(9.8, -6.9, 0.0, 1.0);

    ursa_major[26] = V4(9.8, -6.9, 0.0, 1.0);
    ursa_major[27] = V4(10.2, -8.2, 0.0, 1.0);

    ursa_major[28] = V4(9.8, -6.9, 0.0, 1.0);
    ursa_major[29] = V4(11.0, -7.8, 0.0, 1.0);

    ursa_major[30] = V4(2.1, -1.7, 0.0, 1.0);
    ursa_major[31] = V4(2.6, -4.0, 0.0, 1.0);

    for (int i = 0; i < ursa_major_vert*2; i++){
        ursa_major[i] = Translation(-2, -2, -2.5)*Scaling(.2, .2, .2)*ursa_major[i];
    }


    int orion_vert = 18; 
    vec4* orion = (vec4*) calloc(orion_vert*2, sizeof(vec4));

    orion[0] = V4(2.0, 3.3, 0.0, 1.0);
    orion[1] = V4(1.0, 3.0, 0.0, 1.0);

    orion[2] = V4(1.0, 3.0, 0.0, 1.0);
    orion[3] = V4(0.0, 0.0, 0.0, 1.0);

    orion[4] = V4(0.0, 0.0, 0.0, 1.0);
    orion[5] = V4(0.0, -1.0, 0.0, 1.0);

    orion[6] = V4(0.0, -1.0, 0.0, 1.0);
    orion[7] = V4(1.4, -4.0, 0.0, 1.0);

    orion[8] = V4(1.4, -4.0, 0.0, 1.0);
    orion[9] = V4(2.1, -9.5, 0.0, 1.0);

    orion[10] = V4(2.1, -9.5, 0.0, 1.0);
    orion[11] = V4(1.8, -13.8, 0.0, 1.0);

    orion[12] = V4(1.8, -13.8, 0.0, 1.0);
    orion[13] = V4(6.0, -13.0, 0.0, 1.0);

    orion[14] = V4(6.0, -13.0, 0.0, 1.0);
    orion[15] = V4(4.8, -9.0, 0.0, 1.0);

    orion[16] = V4(4.8, -9.0, 0.0, 1.0);
    orion[17] = V4(3.5, -9.4, 0.0, 1.0);

    orion[18] = V4(3.5, -9.4, 0.0, 1.0);
    orion[19] = V4(2.1, -9.5, 0.0, 1.0);

    orion[20] = V4(4.8, -9.0, 0.0, 1.0);
    orion[21] = V4(5.2, -4.6, 0.0, 1.0);

    orion[22] = V4(5.2, -4.6, 0.0, 1.0);
    orion[23] = V4(4.5, -2.3, 0.0, 1.0);

    orion[24] = V4(4.5, -2.3, 0.0, 1.0);
    orion[25] = V4(1.4, -4.0, 0.0, 1.0);

    orion[26] = V4(5.2, -4.6, 0.0, 1.0);
    orion[27] = V4(11.8, -3.9, 0.0, 1.0);

    orion[28] = V4(10.6, -0.2, 0.0, 1.0);
    orion[29] = V4(11.2, -2.2, 0.0, 1.0);

    orion[30] = V4(11.2, -2.2, 0.0, 1.0);
    orion[31] = V4(11.8, -3.9, 0.0, 1.0);

    orion[32] = V4(11.8, -3.9, 0.0, 1.0);
    orion[33] = V4(11.5, -7.0, 0.0, 1.0);

    orion[34] = V4(11.5, -7.0, 0.0, 1.0);
    orion[35] = V4(10.9, -7.4, 0.0, 1.0);

    mat4 orion_trans = Translation(-4, 3, -2.5)*Scaling(.2, .2, .2);
    for (int i = 0; i < orion_vert*2; i++){
        orion[i] = orion_trans*orion[i];
    }

    const int orion_belt_triangle_mesh_vert_1 = 12;
    vec4 orions_belt_1[orion_belt_triangle_mesh_vert_1] = {
                                                //1, 2, 3, 4 on paper
                                                //2.135 --> 2.065
                                                //used slope intercept formula to get these
                                                V4(2.036, -9.000, 0.0, 1.0), V4(2.100, -9.500, 0.0, 1.0), V4(3.5, -9.4, 0.0, 1.0),
                                                V4(2.065, -10.0, 0.0, 1.0), V4(3.5, -9.4, 0.0, 1.0), V4(2.100, -9.500, 0.0, 1.0),
                                                V4(3.500, -9.400, 0.0, 1.0), V4(4.950, -9.500, 0.0, 1.0), V4(4.800, -9.000, 0.0, 1.0),
                                                V4(3.500, -9.400, 0.0, 1.0), V4(4.800, -9.000, 0.0, 1.0), V4(4.845, -8.500, 0.0, 1.0),
                                                };
    for (int i = 0; i < orion_belt_triangle_mesh_vert_1; i++){
        orions_belt_1[i] = orion_trans*orions_belt_1[i];
    } 

    const int orion_belt_triangle_mesh_vert_2 = 6;
    vec4 orions_belt_2[orion_belt_triangle_mesh_vert_2] = {
                                                //5, 6 on paper
                                                V4(2.036, -9.000, 0.0, 1.0), V4(3.500, -9.400, 0.0, 1.0), V4(4.845, -8.500, 0.0, 1.0),
                                                V4(2.065, -10.0, 0.0, 1.0), V4(3.500, -9.400, 0.0, 1.0), V4(4.950, -9.500, 0.0, 1.0)
                                                };
    for (int i = 0; i < orion_belt_triangle_mesh_vert_2; i++){
        orions_belt_2[i] = orion_trans*orions_belt_2[i];
    } 
    //END CONSTELLATIONS


    //for long hand of clock

    int min_hand_vert = 4;
    vec4* min_hand = (vec4*) calloc(min_hand_vert*2, sizeof(vec4));

    min_hand[0] = V4(0.0, 0.0, 0.0, 1.0);
    min_hand[1] = V4(2.0, 16.0, 0.0, 1.0);

    min_hand[2] = V4(2.0, 16.0, 0.0, 1.0);
    min_hand[3] = V4(0.0, 18.0, 0.0, 1.0);

    min_hand[4] = V4(0.0, 18.0, 0.0, 1.0);
    min_hand[5] = V4(-2.0, 16.0, 0.0, 1.0);

    min_hand[6] = V4(-2.0, 16.0, 0.0, 1.0);
    min_hand[7] = V4(0.0, 0.0, 0.0, 1.0);

    mat4 clock_hand_trans = Translation(0.0, 0.0, -2);

    for (int i = 0; i < min_hand_vert*2; i++){
        min_hand[i] = clock_hand_trans*min_hand[i];
    } 

    
    //start red triangles
    int num_clock_red_tri = 15; 
    vec4* red_tri_min_mesh = (vec4*) calloc(num_clock_red_tri*3, sizeof(vec4));

    for (int tri = 1; tri <= 15; tri += 1){
        int index = (tri-1)*3; 
        red_tri_min_mesh[index] = V4((tri/-8.0), tri, 0.0, 1.0); //bottom left corner
        red_tri_min_mesh[index+1] = V4((tri/8.0), tri, 0.0, 1.0); //bottom right corner
        red_tri_min_mesh[index+2] = V4((tri+1)/-8.0, tri+1, 0.0, 1.0); //top left corner
    }

    for (int i = 0; i < num_clock_red_tri*3; i++){
        red_tri_min_mesh[i] = clock_hand_trans*red_tri_min_mesh[i];
    } 

    //end red triangles

    //start blue triangles
    //same number of blue and red triangles
     int num_clock_blue_tri = 15; 
    vec4* blue_tri_min_mesh = (vec4*) calloc(num_clock_blue_tri*3, sizeof(vec4));

    for (int tri = 1; tri <= 15; tri += 1){
        int index = (tri-1)*3; 
        blue_tri_min_mesh[index] = V4((tri/8.0), tri, 0.0, 1.0); //bottom right corner
        blue_tri_min_mesh[index+1] = V4((tri+1)/8.0, tri+1, 0.0, 1.0); //top right corner
        blue_tri_min_mesh[index+2] = V4((tri+1)/-8.0, tri+1, 0.0, 1.0); //top left corner
    }

    for (int i = 0; i < num_clock_red_tri*3; i++){
        blue_tri_min_mesh[i] = clock_hand_trans*blue_tri_min_mesh[i];
    } 
    //end blue triangles

    //start green triangles
    int num_clock_green_tri = 2;
    vec4* green_tri_min_mesh = (vec4*) calloc(num_clock_green_tri*3, sizeof(vec4));
    //top and bottom of clock hand
    green_tri_min_mesh[0] = clock_hand_trans*V4(0.0, 0.0, 0.0, 1.0);
    green_tri_min_mesh[1] = blue_tri_min_mesh[0];
    green_tri_min_mesh[2] = red_tri_min_mesh[0];

    green_tri_min_mesh[3] = red_tri_min_mesh[14*3+2];
    green_tri_min_mesh[4] = blue_tri_min_mesh[14*3+1];
    green_tri_min_mesh[5] = clock_hand_trans*V4(0.0, 18.0, 0.0, 1.0);
    //end green triangles

    //make hour hand and apply transformations to it
    mat4 hour_hand_trans = Scaling(1.0, .5, 1.0); //make hour hand smaller
    mat4 depth_trans = Translation(0.0, 0.0, -.1);

    vec4* hour_hand = (vec4*) calloc(min_hand_vert*2, sizeof(vec4));
    for (int i = 0; i < min_hand_vert*2; i++){
        hour_hand[i] = depth_trans*hour_hand_trans*min_hand[i];
    } 

    vec4* red_tri_hour_mesh = (vec4*) calloc(num_clock_red_tri*3, sizeof(vec4));
    for (int i = 0; i < num_clock_red_tri*3; i++){
        red_tri_hour_mesh[i] = depth_trans*hour_hand_trans*red_tri_min_mesh[i];
    }

    vec4* blue_tri_hour_mesh = (vec4*) calloc(num_clock_blue_tri*3, sizeof(vec4));
    for (int i = 0; i < num_clock_blue_tri*3; i++){
       blue_tri_hour_mesh[i] = depth_trans*hour_hand_trans*blue_tri_min_mesh[i];
    }

    vec4* green_tri_hour_mesh = (vec4*) calloc(num_clock_green_tri*3, sizeof(vec4));
    for (int i = 0; i < num_clock_green_tri*3; i++){
        green_tri_hour_mesh[i] = depth_trans*hour_hand_trans*green_tri_min_mesh[i];
    }

    int speed = 25;
    bool blinnphong = 0;
    int intensity = 45;

    vec4* rot_red_tri_hour = (vec4*) calloc(num_clock_red_tri*3, sizeof(vec4));
    vec4* rot_blue_tri_hour = (vec4*) calloc(num_clock_blue_tri*3, sizeof(vec4));
    vec4* rot_green_tri_hour = (vec4*) calloc(num_clock_green_tri*3, sizeof(vec4));
    vec4* rot_hour_hand = (vec4*) calloc(min_hand_vert*2, sizeof(vec4));

    vec4* rot_red_tri_min = (vec4*) calloc(num_clock_red_tri*3, sizeof(vec4));
    vec4* rot_blue_tri_min = (vec4*) calloc(num_clock_blue_tri*3, sizeof(vec4));
    vec4* rot_green_tri_min = (vec4*) calloc(num_clock_green_tri*3, sizeof(vec4));
    vec4* rot_min_hand = (vec4*) calloc(min_hand_vert*2, sizeof(vec4));

    bool first = 0;

    //get random stars in background

    int num_stars = 100;
    vec4* background_stars = (vec4*) calloc(100, sizeof(vec4));

    for (int i = 0; i < num_stars; i++){
        int r_x = rand() % 100;
        int r_y = rand() % 100;

        double x = (r_x-50)*(14.0/50);
        double y = (r_y-50)*(14.0/50);

        background_stars[i] = V4(x, y, -5.0, 1.0);
    }

    //BEGIN ANIMATION
    while (begin_frame()) {
        camera_move(&renderer);
            
        { // imgui
            { // reset
                static Camera3D _camera_0 = renderer;
                if (imgui_button("reset", 'r')) {
                    renderer = _camera_0;
                    time = 0;
                    first = false;
                }
            }
            imgui_checkbox("playing", &playing, 'p');
            //imgui_checkbox("shadow_on_ring", &shadow_on_ring, 's');
            imgui_slider("speed", &speed, 5, 100, 'j', 'k');
            imgui_slider("intensity", &intensity, 5, 100, 'u', 'i');

            imgui_checkbox("Blinn-Phong", &blinnphong, 'b');

        }
            
            mat4 P = camera_get_P(&renderer);
            mat4 V = camera_get_V(&renderer);
            
            vec3 x_renderer, y_renderer, z_renderer, o_renderer;
            {
                camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
            }
            
            //draw big_dipper_triangle_mesh
            basic_draw(GL_TRIANGLES, P*V, big_dipper_mesh_vertices, big_dipper_triangle_mesh, monokai.blue);

            //draw big_dipper second; this makes the lines go over the triangles
            basic_draw(GL_LINES,P*V, big_dipper_vert*2, big_dipper, monokai.white);
            
            //draw ursa_major
            basic_draw(GL_LINES, P*V, ursa_major_vert*2, ursa_major, monokai.white);

            //draw orions_belt
            basic_draw(GL_TRIANGLES, P*V, orion_belt_triangle_mesh_vert_1, orions_belt_1, monokai.brown);
            basic_draw(GL_TRIANGLES, P*V, orion_belt_triangle_mesh_vert_2, orions_belt_2, monokai.red);

            //draw orion
            basic_draw(GL_LINES, P*V, orion_vert*2, orion, monokai.white);

            //draw background stars
            basic_draw(GL_POINTS, P*V, num_stars, background_stars, monokai.white, 3.0);

            //draw clock hands over constellations and rotate the clock hands

            //draw hour clock hand first so minute clock hand is in front of hour clock hand
            
            //hour hand should rotate at 1/60th the speed of the minute hand

            mat4 hour_hand_rot; 
            mat4 minute_hand_rot; 

            if (!first){
                hour_hand_rot = Identity4x4; 
                minute_hand_rot = Identity4x4;

                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_red_tri_hour[i] = hour_hand_rot*red_tri_hour_mesh[i];
                } 
                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_blue_tri_hour[i] = hour_hand_rot*blue_tri_hour_mesh[i];
                } 
                for (int i = 0; i < num_clock_green_tri*3; i++){
                    rot_green_tri_hour[i] = hour_hand_rot*green_tri_hour_mesh[i];
                } 
                for (int i = 0; i < min_hand_vert*2; i++){
                    rot_hour_hand[i] = hour_hand_rot*hour_hand[i];
                } 
                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_red_tri_min[i] = minute_hand_rot*red_tri_min_mesh[i];
                } 
                for (int i = 0; i < num_clock_blue_tri*3; i++){
                    rot_blue_tri_min[i] = minute_hand_rot*blue_tri_min_mesh[i];
                } 
                for (int i = 0; i < num_clock_green_tri*3; i++){
                    rot_green_tri_min[i] = minute_hand_rot*green_tri_min_mesh[i];
                } 
                for (int i = 0; i < min_hand_vert*2; i++){
                    rot_min_hand[i] = minute_hand_rot*min_hand[i];
                }
            }

            if (playing){
                //CLAMP RESTRICTS SPEED TO BETWEEN 0 and 50
                first = true;

                // draw hour hand
                hour_hand_rot = RotationZ((-time/(30.0*60.0))*(CLAMP(speed, 0, 50))); 
                

                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_red_tri_hour[i] = hour_hand_rot*red_tri_hour_mesh[i];
                } 
                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_blue_tri_hour[i] = hour_hand_rot*blue_tri_hour_mesh[i];
                } 
                for (int i = 0; i < num_clock_green_tri*3; i++){
                    rot_green_tri_hour[i] = hour_hand_rot*green_tri_hour_mesh[i];
                } 
                
                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_red_tri_hour, monokai.red);
                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_blue_tri_hour, monokai.blue);
                basic_draw(GL_TRIANGLES, P*V, num_clock_green_tri*3, rot_green_tri_hour, monokai.green);

                for (int i = 0; i < min_hand_vert*2; i++){
                    rot_hour_hand[i] = hour_hand_rot*hour_hand[i];
                } 

                //draw minute hand
                minute_hand_rot = RotationZ((-time/30.0)*(CLAMP(speed, 0, 50))); 
                basic_draw(GL_LINES, P*V, min_hand_vert*2, rot_hour_hand, monokai.white);

                for (int i = 0; i < num_clock_red_tri*3; i++){
                    rot_red_tri_min[i] = minute_hand_rot*red_tri_min_mesh[i];
                } 
                for (int i = 0; i < num_clock_blue_tri*3; i++){
                    rot_blue_tri_min[i] = minute_hand_rot*blue_tri_min_mesh[i];
                } 
                for (int i = 0; i < num_clock_green_tri*3; i++){
                    rot_green_tri_min[i] = minute_hand_rot*green_tri_min_mesh[i];
                } 
                
                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_red_tri_min, monokai.red);
                basic_draw(GL_TRIANGLES, P*V, num_clock_blue_tri*3, rot_blue_tri_min, monokai.blue);
                basic_draw(GL_TRIANGLES, P*V, num_clock_green_tri*3, rot_green_tri_min, monokai.green);

                for (int i = 0; i < min_hand_vert*2; i++){
                    rot_min_hand[i] = minute_hand_rot*min_hand[i];
                } 

                basic_draw(GL_LINES, P*V, min_hand_vert*2, rot_min_hand, monokai.white);
            }
            else{
                //do not apply rotation matrices if not playing
                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_red_tri_hour, monokai.red);
                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_blue_tri_hour, monokai.blue);
                basic_draw(GL_TRIANGLES, P*V, num_clock_green_tri*3, rot_green_tri_hour, monokai.green);
                basic_draw(GL_LINES, P*V, min_hand_vert*2, rot_hour_hand, monokai.white);

                basic_draw(GL_TRIANGLES, P*V, num_clock_red_tri*3, rot_red_tri_min, monokai.red);
                basic_draw(GL_TRIANGLES, P*V, num_clock_blue_tri*3, rot_blue_tri_min, monokai.blue);
                basic_draw(GL_TRIANGLES, P*V, num_clock_green_tri*3, rot_green_tri_min, monokai.green);
                basic_draw(GL_LINES, P*V, min_hand_vert*2, rot_min_hand, monokai.white);
            }

            //draw numbers
            double step = 2*3.1415926/24.0;
            vec2 start = V2(planet_positions[0].x, planet_positions[0].y) + V2(0, screen_height_world+1.);

            //24 because of concept art
            char *nums[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                            "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"};
            mat2 rot = {};

            for (int i = 0; i < 24; i++){
                rot = {cos(i*step), sin(i*step), 
                        -sin(i*step), cos(i*step)};

                basic_text(P*V, nums[i], rot*start);
            }

        { // ray marching
            
            
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, num_vertices * 3 * sizeof(double), vertex_positions, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
            glEnableVertexAttribArray(0);

            glUseProgram(shader_program);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_triangles * sizeof(int), triangle_indices, GL_DYNAMIC_DRAW);

            shader_set_uniform_vec2(shader_program, "iResolution", window_get_dimensions_in_pixels());
            {   
                //set time variable in shader
                shader_set_uniform_double(shader_program, "time", time);

                shader_set_uniform_double(shader_program, "renderer_angle_of_view", renderer.angle_of_view);
                vec3 x_renderer, y_renderer, z_renderer, o_renderer; {
                    camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
                }
                
                vec3* planet_change_theta = (vec3*) calloc(1+ num_planets/3, sizeof(vec3));

                //adjust speed 
                planet_change_theta[0] = planet_change_theta_orig[0]*speed;
                planet_change_theta[1] = planet_change_theta_orig[1]*speed;
                planet_change_theta[2] = planet_change_theta_orig[2]*speed;
                planet_change_theta[3] = planet_change_theta_orig[3]*speed;

                //set variables inside of fragment shader
                shader_set_uniform_vec3(shader_program, "x_renderer", x_renderer);
                shader_set_uniform_vec3(shader_program, "y_renderer", y_renderer);
                shader_set_uniform_vec3(shader_program, "z_renderer", z_renderer);
                shader_set_uniform_vec3(shader_program, "o_renderer", o_renderer);
                shader_set_uniform_array_vec3(shader_program, "planet_positions", num_planets, planet_positions);
                shader_set_uniform_array_vec3(shader_program, "planet_colors", num_planets, planet_colors);
                shader_set_uniform_array_vec3(shader_program, "planet_radii", num_planets/3, planet_radii);
                shader_set_uniform_array_vec3(shader_program, "planet_change_theta", num_planets, planet_change_theta);
                shader_set_uniform_array_vec3(shader_program, "asteroid_positions", num_asteroids, asteroid_positions);
                shader_set_uniform_double(shader_program, "asteroid_radius", asteroid_radius);
                shader_set_uniform_int(shader_program, "speed", speed);
                shader_set_uniform_bool(shader_program, "blinnphong", blinnphong);
                shader_set_uniform_bool(shader_program, "shadow_on_ring", shadow_on_ring);
                shader_set_uniform_int(shader_program, "intensity", intensity);


            }

            
            glDrawElements(GL_TRIANGLES, 3 * num_triangles, GL_UNSIGNED_INT, NULL);
        }

        if (playing || input.key_pressed['.']) {
            time += .0167;
        }
    }
}

int main() {
    finalproject();
    return 0;
}

