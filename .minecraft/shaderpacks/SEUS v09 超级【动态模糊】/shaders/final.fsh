#version 120

/*

Settings by Sonic Ether
Bokeh Depth-of-Field by Sonic Ether
God Rays by Blizzard
Bloom shader by CosmicSpore (Modified from original source: http://myheroics.wordpress.com/2008/09/04/glsl-bloom-shader/)
Cross-Processing by Sonic Ether.
High Desaturation effect by Sonic Ether
HDR by Sonic Ether
Glare by Sonic Ether
Shaders 2.0 port of Yourself's Cell Shader, port by an anonymous user.
Bug Fixes by Kool_Kat.

*/




// Place two leading Slashes in front of the following '#define' lines in order to disable an option.
// MOTIONBLUR, HDR, and BOKEH_DOF are very beta shaders. Use at risk of weird results.
// MOTIONBLUR and BOKEH_DOF are not compatable with eachother. Shaders break when you enable both.
// GLARE is still a work in progress.
// BLOOM is currently broken.





//#define BOKEH_DOF
#define GODRAYS
#define GODRAYS_EXPOSURE 0.10
#define GODRAYS_SAMPLES 6
#define GODRAYS_DECAY 0.99
#define GODRAYS_DENSITY 0.10
//#define LENS
//#define LENS_POWER 0.36
//#define GLARE
//#define GLARE_AMOUNT 0.35
//#define GLARE_RANGE 2.0
//#define GLARE2							//second pass of glare shader. More realistic light scattering.
//#define CEL_SHADING
//#define CEL_SHADING_THRESHOLD 0.4
//#define CEL_SHADING_THICKNESS 0.004




//#define VINTAGE
#define VIGNETTE
#define VIGNETTE_STRENGTH 1.30
#define CROSSPROCESS
#define TONEMAP						//Unfinished broken feature
#define BRIGHTMULT 1.1               	// 1.0 = default brightness. Higher values mean brighter. 0 would be black.
#define DARKMULT 0.08						// 0.0 = normal image. Higher values will darken dark colors.
#define COLOR_BOOST	0.00				// 0.0 = normal saturation. Higher values mean more saturated image.
#define MOTIONBLUR
#define MOTIONBLUR_AMOUNT 1.0
#define GAMMA 1.28f							//1.0 is default brightness. lower values will brighten image, higher values will darken image	





// DOF Constants - DO NOT CHANGE
// HYPERFOCAL = (Focal Distance ^ 2)/(Circle of Confusion * F Stop) + Focal Distance
#ifdef USE_DOF
const float HYPERFOCAL = 3.132;
const float PICONSTANT = 3.14159;
#endif





//uniform sampler2D texture;
uniform sampler2D gdepth;
uniform sampler2D composite;
uniform sampler2D gaux1; // red is our motion blur mask. If red == 1, don't blur

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferPreviousProjection;

uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferPreviousModelView;

uniform vec3 cameraPosition;
uniform vec3 previousCameraPosition;

uniform vec3 sunPosition;

uniform float worldTime;
uniform float aspectRatio;
uniform float near;
uniform float far;
uniform float viewWidth;
uniform float viewHeight;
uniform float rainStrength;

varying vec4 texcoord;



//Land/sky mask
float land = texture2D(gaux1, texcoord.st).b;

//Raining
float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;


// Standard depth function.
float getDepth(vec2 coord) {
    return 2.0 * near * far / (far + near - (2.0 * texture2D(gdepth, coord).x - 1.0) * (far - near));
}
float eDepth(vec2 coord) {
	return texture2D(gdepth, coord).x;
}

//Calculate Time of Day

	float timefract = worldTime;

	float TimeSunrise  = ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0) + (1.0 - (clamp(timefract, 0.0, 4000.0)/4000.0));
	float TimeNoon     = ((clamp(timefract, 0.0, 4000.0)) / 4000.0) - ((clamp(timefract, 8000.0, 12000.0) - 8000.0) / 4000.0);
	float TimeSunset   = ((clamp(timefract, 8000.0, 12000.0) - 8000.0) / 4000.0) - ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0);
	float TimeMidnight = ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0) - ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0);



#ifdef BOKEH_DOF

const float blurclamp = 0.014;  // max blur amount
const float bias = 0.3;	//aperture - bigger values for shallower depth of field

#endif

#ifdef GODRAYS

vec3 sunPos = sunPosition;



	float addGodRays(in float nc, in vec2 tx, in float noise, in float noise2, in float noise3, in float noise4, in float noise5, in float noise6, in float noise7, in float noise8, in float noise9) {
			float GDTimeMult = 0.0f;
			if (sunPos.z > 0.0f) {
				sunPos.z = -sunPos.z;
				sunPos.x = -sunPos.x;
				sunPos.y = -sunPos.y;
				GDTimeMult = TimeMidnight;	
			} else {
				GDTimeMult = TimeSunrise + TimeNoon + TimeSunset;
			}
			vec2 lightPos = sunPos.xy / -sunPos.z;
			lightPos.y *= 1.39f;
			lightPos.x *= 0.76f;
			lightPos = (lightPos + 1.0f)/2.0f;
			//vec2 coord = tx;
			vec2 delta = (tx - lightPos) * GODRAYS_DENSITY / float(2.0);
			delta *= -sunPos.z*0.01f;
			//delta *= -sunPos.z*0.01;
			float decay = -sunPos.z / 100.0f;
				 // decay *= -sunPos.z*0.01;
			float colorGD = 0.0f;
			
			for (int i = 0; i < 3; i++) {
			
			if (texcoord.s > 1.0f || texcoord.s < 0.0f || texcoord.t > 1.0f || texcoord.t < 0.0f) {
				break;
			}
			
				
				float sample = 0.0f;

					sample = 1.0f - texture2D(gaux1, tx + vec2(noise*delta.x, noise*delta.y)).g;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise2*delta.x, noise2*delta.y)).g;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise3*delta.x, noise3*delta.y)).g;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise4*delta.x, noise4*delta.y)).g;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise5*delta.x, noise5*delta.y)).g;
					/*
					sample += 1.0 - texture2D(gaux1, tx + vec2(noise6*delta.x, noise6*delta.y)).b;
					sample += 1.0 - texture2D(gaux1, tx + vec2(noise7*delta.x, noise7*delta.y)).b;
					sample += 1.0 - texture2D(gaux1, tx + vec2(noise8*delta.x, noise8*delta.y)).b;
					sample += 1.0 - texture2D(gaux1, tx + vec2(noise9*delta.x, noise9*delta.y)).b;
				*/
				sample *= decay;

					colorGD += sample;
					decay *= GODRAYS_DECAY;
					tx -= delta;
			}
			
			//float bubble = distance(vec2(delta.x*aspectRatio, delta.y), vec2(0.0f, 0.0f))*8.0f;
				 // bubble = clamp(bubble, 0.0f, 1.0f);
				 // bubble = 1.0f - bubble;
				  
			return (nc + GODRAYS_EXPOSURE * (colorGD))*GDTimeMult;
	}
#endif 

#ifdef CEL_SHADING
	float getCellShaderFactor(vec2 coord) {
    float d = getDepth(coord);
    vec3 n = normalize(vec3(getDepth(coord+vec2(CEL_SHADING_THICKNESS,0.0))-d,getDepth(coord+vec2(0.0,CEL_SHADING_THICKNESS))-d , CEL_SHADING_THRESHOLD));
    //clamp(n.z*3.0,0.0,1.0);
    return n.z; 
	}
#endif


// Main --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void main() {

	vec4 color = texture2D(composite, texcoord.st);
	
	
//Common variables

	float depth = eDepth(texcoord.xy);
	const vec2 Texcoord2 = texcoord.st;
	

const float noiseamp = 5.5f;



						const float width3 = 2.0f;
						const float height3 = 2.0f;
						float noiseX3 = ((fract(1.0f-Texcoord2.s*(width3/2.0f))*0.25f)+(fract(Texcoord2.t*(height3/2.0f))*0.75f))*2.0f-1.0f;

						
							noiseX3 = clamp(fract(sin(dot(Texcoord2 ,vec2(18.9898f,28.633f))) * 4378.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX3 *= (0.10f*noiseamp);

						const float width2 = 1.0f;
						const float height2 = 1.0f;
						float noiseX2 = ((fract(1.0f-Texcoord2.s*(width2/2.0f))*0.25f)+(fract(Texcoord2.t*(height2/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY2 = ((fract(1.0f-Texcoord2.s*(width2/2.0f))*0.75f)+(fract(Texcoord2.t*(height2/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX2 = clamp(fract(sin(dot(Texcoord2 ,vec2(12.9898f,78.233f))) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY2 = clamp(fract(sin(dot(Texcoord2 ,vec2(12.9898f,78.233f)*2.0f)) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX2 *= (0.10f*noiseamp);
						noiseY2 *= (0.10f*noiseamp);
						

						const float width4 = 3.0f;
						const float height4 = 3.0f;
						float noiseX4 = ((fract(1.0f-Texcoord2.s*(width4/2.0f))*0.25f)+(fract(Texcoord2.t*(height4/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY4 = ((fract(1.0f-Texcoord2.s*(width4/2.0f))*0.75f)+(fract(Texcoord2.t*(height4/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX4 = clamp(fract(sin(dot(Texcoord2 ,vec2(16.9898f,38.633f))) * 41178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY4 = clamp(fract(sin(dot(Texcoord2 ,vec2(21.9898f,66.233f)*2.0f)) * 9758.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX4 *= (0.10f*noiseamp);
						noiseY4 *= (0.10f*noiseamp);	

						const float width5 = 4.0f;
						const float height5 = 4.0f;
						float noiseX5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.25f)+(fract(Texcoord2.t*(height5/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.75f)+(fract(Texcoord2.t*(height5/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX5 = clamp(fract(sin(dot(Texcoord2 ,vec2(11.9898f,68.633f))) * 21178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY5 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898f,71.233f)*2.0f)) * 6958.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX5 *= (0.10f*noiseamp);
						noiseY5 *= (0.10f*noiseamp);							
						
						const float width6 = 4.0f;
						const float height6 = 4.0f;
						float noiseX6 = ((fract(1.0f-Texcoord2.s*(width6/2.0f))*0.25f)+(fract(Texcoord2.t*(height6/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY6 = ((fract(1.0f-Texcoord2.s*(width6/2.0f))*0.75f)+(fract(Texcoord2.t*(height6/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX6 = clamp(fract(sin(dot(Texcoord2 ,vec2(21.9898f,78.633f))) * 29178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY6 = clamp(fract(sin(dot(Texcoord2 ,vec2(36.9898f,81.233f)*2.0f)) * 16958.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX6 *= (0.10f*noiseamp);
						noiseY6 *= (0.10f*noiseamp);						
						


#ifdef BOKEH_DOF
	
	if (depth > 0.9999f) {
		depth = 1.0f;
	}
	

	float cursorDepth = eDepth(vec2(0.5f, 0.5f));
	
	if (cursorDepth > 0.9999f) {
		cursorDepth = 1.0f;
	}
	
	
	vec2 aspectcorrect = vec2(1.0, aspectRatio) * 1.5;
	
	float factor = (depth - cursorDepth);
	 
	vec2 dofblur = (vec2 (clamp( factor * bias, -blurclamp, blurclamp )))*0.6;

	
	

	vec4 col = vec4(0.0);
	col += texture2D(composite, texcoord.st);
	
	col += texture2D(composite, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur);
	col += texture2D(composite, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur);	
	col += texture2D(composite, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur);
	
	col += texture2D(composite, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(composite, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(composite, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(composite, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(composite, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(composite, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9);	
	col += texture2D(composite, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9);	
	
	col += texture2D(composite, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.7);
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7);
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(composite, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7);
			 
	col += texture2D(composite, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4);
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4);
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(composite, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4);	

	color = col/41;
	

#endif

/*
#ifdef USE_DOF
	float depth = getDepth(texcoord.st);
	    
	float cursorDepth = getDepth(vec2(0.5, 0.5));
    
    // foreground blur = 1/2 background blur. Blur should follow exponential pattern until cursor = hyperfocal -- Cursor before hyperfocal
    // Blur should go from 0 to 1/2 hyperfocal then clear to infinity -- Cursor @ hyperfocal.
    // hyperfocal to inifity is clear though dof extends from 1/2 hyper to hyper -- Cursor beyond hyperfocal
    
    float mixAmount = 0.0;
    
    if (depth < cursorDepth) {
    	mixAmount = clamp(2.0 * ((clamp(cursorDepth, 0.0, HYPERFOCAL) - depth) / (clamp(cursorDepth, 0.0, HYPERFOCAL))), 0.0, 1.0);
	} else if (cursorDepth == HYPERFOCAL) {
		mixAmount = 0.0;
	} else {
		mixAmount =  1.0 - clamp((((cursorDepth * HYPERFOCAL) / (HYPERFOCAL - cursorDepth)) - (depth - cursorDepth)) / ((cursorDepth * HYPERFOCAL) / (HYPERFOCAL - cursorDepth)), 0.0, 1.0);
	}
    
    if (mixAmount != 0.0) {
		color = mix(color, getBlurredColor(), mixAmount);
   	}
#endif
*/






#ifdef MOTIONBLUR

	//float depth = texture2D(gdepth, texcoord.st).x;
	
	float noblur = texture2D(gaux1, texcoord.st).r;

	
		if (depth > 0.9999999f) {
		depth = 1.0f;
		}
		
		float depths = 0.0f;
		const float depthspread = 0.5f;
		float dsx = 0.0f;
		float dsxh = dsx;
		float dsy = 0.0f;
		int dsamples = 0;
	
				for (int i = 0; i < 3; ++i) {
				
						for (int i = 0; i < 3; ++i) {
							
							depths += texture2D(gdepth, texcoord.st + vec2(dsx, dsy)).x;
							dsx += 0.01*depthspread;
							dsamples += 1;
						}
				
					dsy += 0.01*depthspread;
					dsx = dsxh;
					
				}
				
				depths /= dsamples;
				depths = clamp(depths, 0.0f, 0.999999f);
	
	
		if (depth < 1.9999999f) {
		vec4 currentPosition = vec4(texcoord.x * 2.0f - 1.0f, texcoord.y * 2.0f - 1.0f, 2.0f * depths - 1.0f, 1.0f);
	
		vec4 fragposition = gbufferProjectionInverse * currentPosition;
		fragposition = gbufferModelViewInverse * fragposition;
		fragposition /= fragposition.w;
		fragposition.xyz += cameraPosition;
	
		vec4 previousPosition = fragposition;
		previousPosition.xyz -= previousCameraPosition;
		previousPosition = gbufferPreviousModelView * previousPosition;
		previousPosition = gbufferPreviousProjection * previousPosition;
		previousPosition /= previousPosition.w;
	
		vec2 velocity = (currentPosition - previousPosition).st * 0.04f * MOTIONBLUR_AMOUNT;
	
		int samples = 0;
		
		int offsetcount = -2;
		
		
		
		velocity = velocity * (1.0 - noblur);
		
		float edge = distance(texcoord.s, 0.5f);
			  edge = max(edge, distance(texcoord.t, 0.5f));
			  edge *= 2.0f;
			  edge = clamp(edge * 7.0f - 6.0f, 0.0f, 1.0f);
			  edge = 1.0f - edge;
		
		
	
		
		vec2 coord = texcoord.st;
		
		


		for (int i = 0; i < 4; ++i) {
		
		/*
			if (coord.s + velocity.x > 1.0 || coord.t + velocity.y > 1.0 || coord.s + velocity.x < 0.0 || coord.t + velocity.y < 0.0) {
				break;
			}
			*/
			
			coord = coord + (offsetcount*velocity*edge);

			color += texture2D(composite, coord - vec2(noiseX2*velocity.x*edge, noiseX2*velocity.y*edge));
			color += texture2D(composite, coord - vec2(noiseY2*velocity.x*edge, noiseY2*velocity.y*edge));
			color += texture2D(composite, coord - vec2(noiseX4*velocity.x*edge, noiseX4*velocity.y*edge));
			color += texture2D(composite, coord - vec2(noiseY4*velocity.x*edge, noiseY4*velocity.y*edge));
			samples += 4;
			
			offsetcount += 1;
			
			coord = Texcoord2;
		
		}
			color = (color/1.0)/samples;
		}
		

	
#endif



#ifdef GODRAYS

	float GR = addGodRays(0.0f, Texcoord2, noiseX3, noiseX4, noiseY4, noiseX2, noiseY2, noiseX5, noiseY5, noiseX6, noiseY6)/2.0;

	float GRr = 1.0 - texture2D(gaux1, texcoord.st).g;
	
	/*
	float GRs  = 1.0 - texture2D(gaux1, vec2(0.55, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.55, 0.45)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.45)).g;

		  GRs /= 3.0;
	*/
	
	GR = pow(GR, 1.4f)*1.5f;
	
	color.r += GR;
	color.g += GR*0.51f;
	color.b += GR*0.32f;
	
	
	/*
	//Adjust brightness of entire screen based on what the center value of GRs is
	color.r = color.r * (1.0 - (GRs * 0.3));
	color.g = color.g * (1.0 - (GRs * 0.35));
	color.b = color.b * (1.0 - (GRs * 0.5));
	
	color.r = clamp(color.r, 0.0, 1.0);
	color.g = clamp(color.g, 0.0, 1.0);
	color.b = clamp(color.b, 0.0, 1.0);
	
	*/
	
#endif



/*
#ifdef BLOOM
	color = color * 0.8;
	color += addBloom(color, texcoord.st);
#endif
*/


#ifdef GLARE

	color = color * 0.8f;
	
	const float radius = 0.002f*GLARE_RANGE;
	const float radiusv = 0.002f;
	const float bloomintensity = 0.1f*GLARE_AMOUNT;
	
	const float glarex = 2.0f;
	const float glaresub = 0.5f;
	
	float bloomnoise = noiseX2*0.9f;
	

	vec4 clr = vec4(0.0f);
	
	//clr += texture2D(composite, texcoord.st);
	
	//horizontal (70 taps)

	clr +=  max(texture2D(composite, texcoord.st + (vec2(5.0f+bloomnoise,5.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(4.0f+bloomnoise,4.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(3.0f+bloomnoise,3.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(2.0f+bloomnoise,2.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(1.0f+bloomnoise,1.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
		//clr += texture2D(composite, texcoord.st + (vec2(0.0f,0.0f))*radius)*6.0f;
		
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-1.0f+bloomnoise,1.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*5.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-2.0f+bloomnoise,2.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-3.0f+bloomnoise,3.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-4.0f+bloomnoise,4.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-5.0f+bloomnoise,5.0+bloomnoise))*radius)*glarex - glaresub, 0.0f)*1.0f;

	//vertical

	clr +=  max(texture2D(composite, texcoord.st + (vec2(5.0+bloomnoise,-5.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(4.0+bloomnoise,-4.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(3.0+bloomnoise,-3.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(2.0+bloomnoise,-2.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(1.0+bloomnoise,-1.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-5.0+bloomnoise,-5.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-4.0+bloomnoise,-4.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-3.0+bloomnoise,-3.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-2.0+bloomnoise,-2.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-1.0+bloomnoise,-1.0f+bloomnoise))*radius)*glarex - glaresub, 0.0f)*5.0f;
	

	
	clr = (clr/10.0f)/2.0f;
	
	const float clrboost = 0.0;
	
	clr.r = clr.r + (clr.r*(clrboost*2.0)) - (clr.g * clrboost) - (clr.b * clrboost);
	clr.g = clr.g + (clr.g*(clrboost*2.0)) - (clr.r * clrboost) - (clr.b * clrboost);
	clr.b = clr.b + (clr.b*(clrboost*2.0)) - (clr.r * clrboost) - (clr.g * clrboost);

	
	color.r = color.r + (clr.r*1.0f)*bloomintensity;
	color.g = color.g + (clr.g*1.0f)*bloomintensity;
	color.b = color.b + (clr.b*1.0f)*bloomintensity;
	color = max(color, 0.0f);
	

#endif




#ifdef GLARE2

	color = color * 0.8f;
	
	const float radius2 = 0.006f*GLARE_RANGE;
	const float radius2v = 0.002f;
	const float bloomintensity2 = 0.08f*GLARE_AMOUNT;
	
	const float glarex2 = 2.0f;
	const float glaresub2 = 0.5f;
	
	float bloomnoise2 = noiseX4*0.9f;	

	vec4 clr2 = vec4(0.0f);
	
	//clr2 += texture2D(composite, texcoord.st);
	
	//horizontal (70 taps)

	clr2 += max(texture2D(composite, texcoord.st + (vec2(5.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(4.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(3.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(2.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(1.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
		//clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0f,0.0f))*radius2)*6.0f;
		
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-1.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-2.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-3.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-4.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-5.0f+bloomnoise2,0.0+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;

	//vertical

	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-5.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-4.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-3.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-2.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-1.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,5.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,4.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,3.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,2.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,1.0f+bloomnoise2))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	

	
	clr2 = (clr2/10.0f)/2.0f;
	
	const float clr2boost = 0.0;
	
	clr2.r = clr2.r + (clr2.r*(clr2boost*2.0)) - (clr2.g * clr2boost) - (clr2.b * clr2boost);
	clr2.g = clr2.g + (clr2.g*(clr2boost*2.0)) - (clr2.r * clr2boost) - (clr2.b * clr2boost);
	clr2.b = clr2.b + (clr2.b*(clr2boost*2.0)) - (clr2.r * clr2boost) - (clr2.g * clr2boost);

	
	color.r = color.r + (clr2.r*1.0f)*bloomintensity2;
	color.g = color.g + (clr2.g*1.0f)*bloomintensity2;
	color.b = color.b + (clr2.b*1.5f)*bloomintensity2;
	color = max(color, 0.0f);
	

#endif





#ifdef VIGNETTE

float dv = distance(texcoord.st, vec2(0.5f, 0.5f));

dv *= VIGNETTE_STRENGTH;

dv = 1.0f - dv;

dv = pow(dv, 0.2f);

dv *= 1.9f;
dv -= 0.9f;

color.r = color.r * dv;
color.g = color.g * dv;
color.b = color.b * dv;

#endif






#ifdef LENS

vec3 sP = sunPosition;

			vec2 lPos = sP.xy / -sP.z;
			lPos.y *= 1.39f;
			lPos.x *= 0.76f;
			lPos = (lPos + 1.0f)/2.0f;
			
			float sunmask = 0.0f;
			float sunstep = -4.5f;
			float masksize = 0.004f;
					
					for (int i = 0; i < 1; ++i) {
					
					if (lPos.s < 0.0 || lPos.s > 1.0 || lPos.t < 0.0 || lPos.t > 1.0) {
						break;
					}
					
					
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(-5.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(-4.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(-3.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(-2.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(-0.0f*masksize, sunstep*masksize)).b;
					sunmask += 1.0f - texture2D(gaux1, lPos).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(1.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(2.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(3.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(4.0f*masksize, sunstep*masksize)).b;
					//sunmask += 1.0f - texture2D(gaux1, lPos + vec2(5.0f*masksize, sunstep*masksize)).b;
					
					sunstep += 1.0f;
					
					}
					sunmask /= 1.0f;
					sunmask *= LENS_POWER;
					sunmask *= 1.0 - rainx;
			
			//Detect if sun is on edge of screen
				float edgemaskx = clamp(distance(lPos.x, 0.5f)*8.0f - 3.0f, 0.0f, 1.0f);
				float edgemasky = clamp(distance(lPos.y, 0.5f)*8.0f - 3.0f, 0.0f, 1.0f);
			
						
						
			////Darken colors if the sun is visible
				float centermask = 1.0 - clamp(distance(lPos.xy, vec2(0.5f, 0.5f))*2.0, 0.0, 1.0);
						centermask = pow(centermask, 1.0f);
						centermask *= sunmask;
			
				color.r *= (1.0 - centermask);
				color.g *= (1.0 - centermask);
				color.b *= (1.0 - centermask);
			
			
			//Adjust global flare settings
				const float flaremultR = 1.0f;
				const float flaremultG = 1.0f;
				const float flaremultB = 1.0f;
			
				float flarescale = 1.0f;
				const float flarescaleconst = 1.0f;
			
			
			//Flare gets bigger at center of screen
			
				flarescale *= (1.0 - centermask);
			

			//Center white flare
			const vec2 flare1scale = vec2(1.7f*flarescale, 1.7f*flarescale);
			const float flare1pow = 12.0f;
			vec2 flare1pos = vec2(lPos.x*aspectRatio*flare1scale.x, lPos.y*flare1scale.y);
			
			
			float flare1 = distance(flare1pos, vec2(texcoord.s*aspectRatio*flare1scale.x, texcoord.t*flare1scale.y));
				  flare1 = 0.5 - flare1;
				  flare1 = clamp(flare1, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  flare1 *= sunmask;
				  flare1 = pow(flare1, 1.8f);
				  
				  flare1 *= flare1pow;
				  
				  	color.r += flare1*0.7f*flaremultR;
					color.g += flare1*0.4f*flaremultG;
					color.b += flare1*0.2f*flaremultB;	
				  			
							
							
			//Center white flare
			const vec2 flare1Bscale = vec2(0.5f*flarescale, 0.5f*flarescale);
			const float flare1Bpow = 6.0f;
			vec2 flare1Bpos = vec2(lPos.x*aspectRatio*flare1Bscale.x, lPos.y*flare1Bscale.y);
			
			
			float flare1B = distance(flare1Bpos, vec2(texcoord.s*aspectRatio*flare1Bscale.x, texcoord.t*flare1Bscale.y));
				  flare1B = 0.5 - flare1B;
				  flare1B = clamp(flare1B, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  flare1B *= sunmask;
				  flare1B = pow(flare1B, 1.8f);
				  
				  flare1B *= flare1Bpow;
				  
				  	color.r += flare1B*0.7f*flaremultR;
					color.g += flare1B*0.2f*flaremultG;
					color.b += flare1B*0.0f*flaremultB;	
				  
				  
			//Wide red flare
			vec2 flare2pos = vec2(lPos.x*aspectRatio*0.2, lPos.y);
			
			float flare2 = distance(flare2pos, vec2(texcoord.s*aspectRatio*0.2, texcoord.t));
				  flare2 = 0.3 - flare2;
				  flare2 = clamp(flare2, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  flare2 *= sunmask;
				  flare2 = pow(flare2, 1.8f);
				  	
					color.r += flare2*1.8f*flaremultR;
					color.g += flare2*0.6f*flaremultG;
					color.b += flare2*0.0f*flaremultB;
					
					
					
			//Wide red flare
			vec2 flare2posB = vec2(lPos.x*aspectRatio*0.2, lPos.y*4.0);
			
			float flare2B = distance(flare2posB, vec2(texcoord.s*aspectRatio*0.2, texcoord.t*4.0));
				  flare2B = 0.3 - flare2B;
				  flare2B = clamp(flare2B, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  flare2B *= sunmask;
				  flare2B = pow(flare2B, 1.8f);
				  	
					color.r += flare2B*1.2f*flaremultR;
					color.g += flare2B*0.5f*flaremultG;
					color.b += flare2B*0.0f*flaremultB;
					
					
					
			//Far blue flare MAIN
			const vec2 flare3scale = vec2(2.0f*flarescale, 2.0f*flarescale);
			const float flare3pow = 0.7f;
			const float flare3fill = 10.0f;
			const float flare3offset = -0.5f;
			vec2 flare3pos = vec2(  ((1.0 - lPos.x)*(flare3offset + 1.0) - (flare3offset*0.5))  *aspectRatio*flare3scale.x,  ((1.0 - lPos.y)*(flare3offset + 1.0) - (flare3offset*0.5))  *flare3scale.y);
			
			
			float flare3 = distance(flare3pos, vec2(texcoord.s*aspectRatio*flare3scale.x, texcoord.t*flare3scale.y));
				  flare3 = 0.5 - flare3;
				  flare3 = clamp(flare3*flare3fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare3 = sin(flare3*1.57075);
				  flare3 *= sunmask;
				  flare3 = pow(flare3, 1.1f);
				  
				  flare3 *= flare3pow;			
				  
				  
				  //subtract from blue flare
				const vec2 flare3Bscale = vec2(1.4f*flarescale, 1.4f*flarescale);
				const float flare3Bpow = 1.0f;
				const float flare3Bfill = 2.0f;
				const float flare3Boffset = -0.65f;
				vec2 flare3Bpos = vec2(  ((1.0 - lPos.x)*(flare3Boffset + 1.0) - (flare3Boffset*0.5))  *aspectRatio*flare3Bscale.x,  ((1.0 - lPos.y)*(flare3Boffset + 1.0) - (flare3Boffset*0.5))  *flare3Bscale.y);
			
			
				float flare3B = distance(flare3Bpos, vec2(texcoord.s*aspectRatio*flare3Bscale.x, texcoord.t*flare3Bscale.y));
					flare3B = 0.5 - flare3B;
					flare3B = clamp(flare3B*flare3Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
					flare3B = sin(flare3B*1.57075);
					flare3B *= sunmask;
					flare3B = pow(flare3B, 0.9f);
				  
					flare3B *= flare3Bpow;
				  
				flare3 = clamp(flare3 - flare3B, 0.0, 10.0);
				  
				  
				  	color.r += flare3*0.0f*flaremultR;
					color.g += flare3*0.3f*flaremultG;
					color.b += flare3*1.0f*flaremultB;

					
					
					
			//Far blue flare MAIN 2
			const vec2 flare3Cscale = vec2(3.2f*flarescale, 3.2f*flarescale);
			const float flare3Cpow = 1.4f;
			const float flare3Cfill = 10.0f;
			const float flare3Coffset = -0.0f;
			vec2 flare3Cpos = vec2(  ((1.0 - lPos.x)*(flare3Coffset + 1.0) - (flare3Coffset*0.5))  *aspectRatio*flare3Cscale.x,  ((1.0 - lPos.y)*(flare3Coffset + 1.0) - (flare3Coffset*0.5))  *flare3Cscale.y);
			
			
			float flare3C = distance(flare3Cpos, vec2(texcoord.s*aspectRatio*flare3Cscale.x, texcoord.t*flare3Cscale.y));
				  flare3C = 0.5 - flare3C;
				  flare3C = clamp(flare3C*flare3Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare3C = sin(flare3C*1.57075);
				  
				  flare3C = pow(flare3C, 1.1f);
				  
				  flare3C *= flare3Cpow;			
				  
				  
				  //subtract from blue flare
				const vec2 flare3Dscale = vec2(2.1f*flarescale, 2.1f*flarescale);
				const float flare3Dpow = 2.7f;
				const float flare3Dfill = 1.4f;
				const float flare3Doffset = -0.05f;
				vec2 flare3Dpos = vec2(  ((1.0 - lPos.x)*(flare3Doffset + 1.0) - (flare3Doffset*0.5))  *aspectRatio*flare3Dscale.x,  ((1.0 - lPos.y)*(flare3Doffset + 1.0) - (flare3Doffset*0.5))  *flare3Dscale.y);
			
			
				float flare3D = distance(flare3Dpos, vec2(texcoord.s*aspectRatio*flare3Dscale.x, texcoord.t*flare3Dscale.y));
					flare3D = 0.5 - flare3D;
					flare3D = clamp(flare3D*flare3Dfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
					flare3D = sin(flare3D*1.57075);
					flare3D = pow(flare3D, 0.9f);
				  
					flare3D *= flare3Dpow;
				  
				flare3C = clamp(flare3C - flare3D, 0.0, 10.0);
				flare3C *= sunmask;
				  
				  	color.r += flare3C*0.4f*flaremultR;
					color.g += flare3C*0.7f*flaremultG;
					color.b += flare3C*1.0f*flaremultB;							
					
					
					
					
					
					
					
					
					
			//far small pink flare
			const vec2 flare4scale = vec2(4.5f*flarescale, 4.5f*flarescale);
			const float flare4pow = 0.3f;
			const float flare4fill = 3.0f;
			const float flare4offset = -0.1f;
			vec2 flare4pos = vec2(  ((1.0 - lPos.x)*(flare4offset + 1.0) - (flare4offset*0.5))  *aspectRatio*flare4scale.x,  ((1.0 - lPos.y)*(flare4offset + 1.0) - (flare4offset*0.5))  *flare4scale.y);
			
			
			float flare4 = distance(flare4pos, vec2(texcoord.s*aspectRatio*flare4scale.x, texcoord.t*flare4scale.y));
				  flare4 = 0.5 - flare4;
				  flare4 = clamp(flare4*flare4fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4 = sin(flare4*1.57075);
				  flare4 *= sunmask;
				  flare4 = pow(flare4, 1.1f);
				  
				  flare4 *= flare4pow;
				  
				  	color.r += flare4*0.6f*flaremultR;
					color.g += flare4*0.0f*flaremultG;
					color.b += flare4*0.8f*flaremultB;							
					
					
					
			//far small pink flare2
			const vec2 flare4Bscale = vec2(7.5f*flarescale, 7.5f*flarescale);
			const float flare4Bpow = 0.4f;
			const float flare4Bfill = 2.0f;
			const float flare4Boffset = 0.0f;
			vec2 flare4Bpos = vec2(  ((1.0 - lPos.x)*(flare4Boffset + 1.0) - (flare4Boffset*0.5))  *aspectRatio*flare4Bscale.x,  ((1.0 - lPos.y)*(flare4Boffset + 1.0) - (flare4Boffset*0.5))  *flare4Bscale.y);
			
			
			float flare4B = distance(flare4Bpos, vec2(texcoord.s*aspectRatio*flare4Bscale.x, texcoord.t*flare4Bscale.y));
				  flare4B = 0.5 - flare4B;
				  flare4B = clamp(flare4B*flare4Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4B = sin(flare4B*1.57075);
				  flare4B *= sunmask;
				  flare4B = pow(flare4B, 1.1f);
				  
				  flare4B *= flare4Bpow;
				  
				  	color.r += flare4B*0.4f*flaremultR;
					color.g += flare4B*0.0f*flaremultG;
					color.b += flare4B*0.8f*flaremultB;						
					
					
					
			//far small pink flare3
			const vec2 flare4Cscale = vec2(37.5f*flarescale, 37.5f*flarescale);
			const float flare4Cpow = 2.0f;
			const float flare4Cfill = 2.0f;
			const float flare4Coffset = -0.3f;
			vec2 flare4Cpos = vec2(  ((1.0 - lPos.x)*(flare4Coffset + 1.0) - (flare4Coffset*0.5))  *aspectRatio*flare4Cscale.x,  ((1.0 - lPos.y)*(flare4Coffset + 1.0) - (flare4Coffset*0.5))  *flare4Cscale.y);
			
			
			float flare4C = distance(flare4Cpos, vec2(texcoord.s*aspectRatio*flare4Cscale.x, texcoord.t*flare4Cscale.y));
				  flare4C = 0.5 - flare4C;
				  flare4C = clamp(flare4C*flare4Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4C = sin(flare4C*1.57075);
				  flare4C *= sunmask;
				  flare4C = pow(flare4C, 1.1f);
				  
				  flare4C *= flare4Cpow;
				  
				  	color.r += flare4C*0.2f*flaremultR;
					color.g += flare4C*0.6f*flaremultG;
					color.b += flare4C*0.8f*flaremultB;						
					
					
					
			//far small pink flare4
			const vec2 flare4Dscale = vec2(67.5f*flarescale, 67.5f*flarescale);
			const float flare4Dpow = 1.0f;
			const float flare4Dfill = 2.0f;
			const float flare4Doffset = -0.35f;
			vec2 flare4Dpos = vec2(  ((1.0 - lPos.x)*(flare4Doffset + 1.0) - (flare4Doffset*0.5))  *aspectRatio*flare4Dscale.x,  ((1.0 - lPos.y)*(flare4Doffset + 1.0) - (flare4Doffset*0.5))  *flare4Dscale.y);
			
			
			float flare4D = distance(flare4Dpos, vec2(texcoord.s*aspectRatio*flare4Dscale.x, texcoord.t*flare4Dscale.y));
				  flare4D = 0.5 - flare4D;
				  flare4D = clamp(flare4D*flare4Dfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4D = sin(flare4D*1.57075);
				  flare4D *= sunmask;
				  flare4D = pow(flare4D, 1.1f);
				  
				  flare4D *= flare4Dpow;
				  
				  	color.r += flare4D*0.2f*flaremultR;
					color.g += flare4D*0.2f*flaremultG;
					color.b += flare4D*0.8f*flaremultB;						
					
					
								
			//far small pink flare5
			const vec2 flare4Escale = vec2(60.5f*flarescale, 60.5f*flarescale);
			const float flare4Epow = 1.0f;
			const float flare4Efill = 3.0f;
			const float flare4Eoffset = -0.3393f;
			vec2 flare4Epos = vec2(  ((1.0 - lPos.x)*(flare4Eoffset + 1.0) - (flare4Eoffset*0.5))  *aspectRatio*flare4Escale.x,  ((1.0 - lPos.y)*(flare4Eoffset + 1.0) - (flare4Eoffset*0.5))  *flare4Escale.y);
			
			
			float flare4E = distance(flare4Epos, vec2(texcoord.s*aspectRatio*flare4Escale.x, texcoord.t*flare4Escale.y));
				  flare4E = 0.5 - flare4E;
				  flare4E = clamp(flare4E*flare4Efill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4E = sin(flare4E*1.57075);
				  flare4E *= sunmask;
				  flare4E = pow(flare4E, 1.1f);
				  
				  flare4E *= flare4Epow;
				  
				  	color.r += flare4E*0.2f*flaremultR;
					color.g += flare4E*0.2f*flaremultG;
					color.b += flare4E*0.6f*flaremultB;					
					
								
								
			//far small pink flare5
			const vec2 flare4Fscale = vec2(20.5f*flarescale, 20.5f*flarescale);
			const float flare4Fpow = 3.0f;
			const float flare4Ffill = 3.0f;
			const float flare4Foffset = -0.4713f;
			vec2 flare4Fpos = vec2(  ((1.0 - lPos.x)*(flare4Foffset + 1.0) - (flare4Foffset*0.5))  *aspectRatio*flare4Fscale.x,  ((1.0 - lPos.y)*(flare4Foffset + 1.0) - (flare4Foffset*0.5))  *flare4Fscale.y);
			
			
			float flare4F = distance(flare4Fpos, vec2(texcoord.s*aspectRatio*flare4Fscale.x, texcoord.t*flare4Fscale.y));
				  flare4F = 0.5 - flare4F;
				  flare4F = clamp(flare4F*flare4Ffill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare4F = sin(flare4F*1.57075);
				  flare4F *= sunmask;
				  flare4F = pow(flare4F, 1.1f);
				  
				  flare4F *= flare4Fpow;
				  
				  	color.r += flare4F*0.6f*flaremultR;
					color.g += flare4F*0.1f*flaremultG;
					color.b += flare4F*0.1f*flaremultB;						
					
					
					
					
					
					
					
					
					
					
					
					
			//
			const vec2 flare5scale = vec2(3.2f*flarescaleconst, 3.2f*flarescaleconst);
			const float flare5pow = 13.4f;
			const float flare5fill = 1.0f;
			const float flare5offset = -2.0f;
			vec2 flare5pos = vec2(  ((1.0 - lPos.x)*(flare5offset + 1.0) - (flare5offset*0.5))  *aspectRatio*flare5scale.x,  ((1.0 - lPos.y)*(flare5offset + 1.0) - (flare5offset*0.5))  *flare5scale.y);
			
			
			float flare5 = distance(flare5pos, vec2(texcoord.s*aspectRatio*flare5scale.x, texcoord.t*flare5scale.y));
				  flare5 = 0.5 - flare5;
				  flare5 = clamp(flare5*flare5fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare5 *= sunmask;
				  flare5 = pow(flare5, 1.9f);
				  
				  flare5 *= flare5pow;
				  
				  	color.r += flare5*0.9f*flaremultR;
					color.g += flare5*0.4f*flaremultG;
					color.b += flare5*0.3f*flaremultB;						
					
				/*	
					
			//Soft blue strip 
			const vec2 flare5Bscale = vec2(0.5f*flarescaleconst, 3.5f*flarescaleconst);
			const float flare5Bpow = 1.4f;
			const float flare5Bfill = 2.0f;
			const float flare5Boffset = -1.9f;
			vec2 flare5Bpos = vec2(  ((1.0 - lPos.x)*(flare5Boffset + 1.0) - (flare5Boffset*0.5))  *aspectRatio*flare5Bscale.x,  ((1.0 - lPos.y)*(flare5Boffset + 1.0) - (flare5Boffset*0.5))  *flare5Bscale.y);
			
			
			float flare5B = distance(flare5Bpos, vec2(texcoord.s*aspectRatio*flare5Bscale.x, texcoord.t*flare5Bscale.y));
				  flare5B = 0.5 - flare5B;
				  flare5B = clamp(flare5B*flare5Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare5B *= sunmask;
				  flare5B = pow(flare5B, 1.4f);
				  
				  flare5B *= flare5Bpow;
				  
				  	color.r += flare5B*0.9f*flaremultR;
					color.g += flare5B*0.3f*flaremultG;
					color.b += flare5B*0.0f*flaremultB;						
					
					*/
					
					
			//close ring flare red
			const vec2 flare6scale = vec2(1.2f*flarescale, 1.2f*flarescale);
			const float flare6pow = 0.2f;
			const float flare6fill = 5.0f;
			const float flare6offset = -1.9f;
			vec2 flare6pos = vec2(  ((1.0 - lPos.x)*(flare6offset + 1.0) - (flare6offset*0.5))  *aspectRatio*flare6scale.x,  ((1.0 - lPos.y)*(flare6offset + 1.0) - (flare6offset*0.5))  *flare6scale.y);
			
			
			float flare6 = distance(flare6pos, vec2(texcoord.s*aspectRatio*flare6scale.x, texcoord.t*flare6scale.y));
				  flare6 = 0.5 - flare6;
				  flare6 = clamp(flare6*flare6fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare6 = pow(flare6, 1.6f);
				  flare6 = sin(flare6*3.1415);
				  flare6 *= sunmask;

				  
				  flare6 *= flare6pow;
				  
				  	color.r += flare6*0.6f*flaremultR;
					color.g += flare6*0.0f*flaremultG;
					color.b += flare6*0.0f*flaremultB;						
					
					
					
			//close ring flare green
			const vec2 flare6Bscale = vec2(1.1f*flarescale, 1.1f*flarescale);
			const float flare6Bpow = 0.2f;
			const float flare6Bfill = 5.0f;
			const float flare6Boffset = -1.9f;
			vec2 flare6Bpos = vec2(  ((1.0 - lPos.x)*(flare6Boffset + 1.0) - (flare6Boffset*0.5))  *aspectRatio*flare6Bscale.x,  ((1.0 - lPos.y)*(flare6Boffset + 1.0) - (flare6Boffset*0.5))  *flare6Bscale.y);
			
			
			float flare6B = distance(flare6Bpos, vec2(texcoord.s*aspectRatio*flare6Bscale.x, texcoord.t*flare6Bscale.y));
				  flare6B = 0.5 - flare6B;
				  flare6B = clamp(flare6B*flare6Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare6B = pow(flare6B, 1.6f);
				  flare6B = sin(flare6B*3.1415);
				  flare6B *= sunmask;

				  
				  flare6B *= flare6Bpow;
				  
				  	color.r += flare6B*0.0f*flaremultR;
					color.g += flare6B*0.4f*flaremultG;
					color.b += flare6B*0.0f*flaremultB;						
					
					
			
			//close ring flare blue
			const vec2 flare6Cscale = vec2(0.9f*flarescale, 0.9f*flarescale);
			const float flare6Cpow = 0.2f;
			const float flare6Cfill = 5.0f;
			const float flare6Coffset = -1.9f;
			vec2 flare6Cpos = vec2(  ((1.0 - lPos.x)*(flare6Coffset + 1.0) - (flare6Coffset*0.5))  *aspectRatio*flare6Cscale.x,  ((1.0 - lPos.y)*(flare6Coffset + 1.0) - (flare6Coffset*0.5))  *flare6Cscale.y);
			
			
			float flare6C = distance(flare6Cpos, vec2(texcoord.s*aspectRatio*flare6Cscale.x, texcoord.t*flare6Cscale.y));
				  flare6C = 0.5 - flare6C;
				  flare6C = clamp(flare6C*flare6Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare6C = pow(flare6C, 1.8f);
				  flare6C = sin(flare6C*3.1415);
				  flare6C *= sunmask;

				  
				  flare6C *= flare6Cpow;
				  
				  	color.r += flare6C*0.0f*flaremultR;
					color.g += flare6C*0.0f*flaremultG;
					color.b += flare6C*0.4f*flaremultB;						
					
					
					
					
			//far red ring

			const vec2 flare7scale = vec2(0.4f*flarescale, 0.4f*flarescale);
			const float flare7pow = 0.2f;
			const float flare7fill = 10.0f;
			const float flare7offset = 2.6f;
			vec2 flare7pos = vec2(  ((1.0 - lPos.x)*(flare7offset + 1.0) - (flare7offset*0.5))  *aspectRatio*flare7scale.x,  ((1.0 - lPos.y)*(flare7offset + 1.0) - (flare7offset*0.5))  *flare7scale.y);
			
			
			float flare7 = distance(flare7pos, vec2(texcoord.s*aspectRatio*flare7scale.x, texcoord.t*flare7scale.y));
				  flare7 = 0.5 - flare7;
				  flare7 = clamp(flare7*flare7fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare7 = pow(flare7, 1.9f);
				  flare7 = sin(flare7*3.1415);
				  flare7 *= sunmask;

				  
				  flare7 *= flare7pow;
				  
				  	color.r += flare7*1.0f*flaremultR;
					color.g += flare7*0.0f*flaremultG;
					color.b += flare7*0.0f*flaremultB;				
					
					
					
			//far blue ring

			const vec2 flare7Dscale = vec2(0.39f*flarescale, 0.39f*flarescale);
			const float flare7Dpow = 0.1f;
			const float flare7Dfill = 10.0f;
			const float flare7Doffset = 2.6f;
			vec2 flare7Dpos = vec2(  ((1.0 - lPos.x)*(flare7Doffset + 1.0) - (flare7Doffset*0.5))  *aspectRatio*flare7Dscale.x,  ((1.0 - lPos.y)*(flare7Doffset + 1.0) - (flare7Doffset*0.5))  *flare7Dscale.y);
			
			
			float flare7D = distance(flare7Dpos, vec2(texcoord.s*aspectRatio*flare7Dscale.x, texcoord.t*flare7Dscale.y));
				  flare7D = 0.5 - flare7D;
				  flare7D = clamp(flare7D*flare7Dfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare7D = pow(flare7D, 1.9f);
				  flare7D = sin(flare7D*3.1415);
				  flare7D *= sunmask;

				  
				  flare7D *= flare7Dpow;
				  
				  	color.r += flare7D*0.0f*flaremultR;
					color.g += flare7D*0.6f*flaremultG;
					color.b += flare7D*0.0f*flaremultB;				
					
					
					
			//far red glow

			const vec2 flare7Bscale = vec2(0.2f*flarescale, 0.2f*flarescale);
			const float flare7Bpow = 0.1f;
			const float flare7Bfill = 2.0f;
			const float flare7Boffset = 2.9f;
			vec2 flare7Bpos = vec2(  ((1.0 - lPos.x)*(flare7Boffset + 1.0) - (flare7Boffset*0.5))  *aspectRatio*flare7Bscale.x,  ((1.0 - lPos.y)*(flare7Boffset + 1.0) - (flare7Boffset*0.5))  *flare7Bscale.y);
			
			
			float flare7B = distance(flare7Bpos, vec2(texcoord.s*aspectRatio*flare7Bscale.x, texcoord.t*flare7Bscale.y));
				  flare7B = 0.5 - flare7B;
				  flare7B = clamp(flare7B*flare7Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare7B = pow(flare7B, 1.9f);
				  flare7B = sin(flare7B*3.1415*0.5);
				  flare7B *= sunmask;

				  
				  flare7B *= flare7Bpow;
				  
				  	color.r += flare7B*1.0f*flaremultR;
					color.g += flare7B*0.0f*flaremultG;
					color.b += flare7B*0.0f*flaremultB;	
			
			
			
			//Edge blue strip 1
			const vec2 flare8scale = vec2(0.3f*flarescale, 40.5f*flarescale);
			const float flare8pow = 0.5f;
			const float flare8fill = 12.0f;
			const float flare8offset = 1.0f;
			vec2 flare8pos = vec2(  ((1.0 - lPos.x)*(flare8offset + 1.0) - (flare8offset*0.5))  *aspectRatio*flare8scale.x,  ((lPos.y)*(flare8offset + 1.0) - (flare8offset*0.5))  *flare8scale.y);
			
			
			float flare8 = distance(flare8pos, vec2(texcoord.s*aspectRatio*flare8scale.x, texcoord.t*flare8scale.y));
				  flare8 = 0.5 - flare8;
				  flare8 = clamp(flare8*flare8fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare8 *= sunmask;
				  flare8 = pow(flare8, 1.4f);
				  
				  flare8 *= flare8pow;
				  flare8 *= edgemaskx;
				  
				  	color.r += flare8*0.0f*flaremultR;
					color.g += flare8*0.3f*flaremultG;
					color.b += flare8*0.8f*flaremultB;					
			
		
		
			//Edge blue strip 1
			const vec2 flare9scale = vec2(0.2f*flarescale, 5.5f*flarescale);
			const float flare9pow = 1.9f;
			const float flare9fill = 2.0f;
			const vec2 flare9offset = vec2(1.0f, 0.0f);
			vec2 flare9pos = vec2(  ((1.0 - lPos.x)*(flare9offset.x + 1.0) - (flare9offset.x*0.5))  *aspectRatio*flare9scale.x,  ((1.0 - lPos.y)*(flare9offset.y + 1.0) - (flare9offset.y*0.5))  *flare9scale.y);
			
			
			float flare9 = distance(flare9pos, vec2(texcoord.s*aspectRatio*flare9scale.x, texcoord.t*flare9scale.y));
				  flare9 = 0.5 - flare9;
				  flare9 = clamp(flare9*flare9fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare9 *= sunmask;
				  flare9 = pow(flare9, 1.4f);
				  
				  flare9 *= flare9pow;
				  flare9 *= edgemaskx;
				  
				  	color.r += flare9*0.2f*flaremultR;
					color.g += flare9*0.4f*flaremultG;
					color.b += flare9*0.9f*flaremultB;		
					
					
					
		//SMALL SWEEPS		///////////////////////////////						
					
					
			//mid orange sweep
			const vec2 flare10scale = vec2(6.0f*flarescale, 6.0f*flarescale);
			const float flare10pow = 1.9f;
			const float flare10fill = 1.1f;
			const float flare10offset = -0.7f;
			vec2 flare10pos = vec2(  ((1.0 - lPos.x)*(flare10offset + 1.0) - (flare10offset*0.5))  *aspectRatio*flare10scale.x,  ((1.0 - lPos.y)*(flare10offset + 1.0) - (flare10offset*0.5))  *flare10scale.y);
			
			
			float flare10 = distance(flare10pos, vec2(texcoord.s*aspectRatio*flare10scale.x, texcoord.t*flare10scale.y));
				  flare10 = 0.5 - flare10;
				  flare10 = clamp(flare10*flare10fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare10 = sin(flare10*1.57075);
				  flare10 *= sunmask;
				  flare10 = pow(flare10, 1.1f);
				  
				  flare10 *= flare10pow;			
				  
				  
				  //subtract
				const vec2 flare10Bscale = vec2(5.1f*flarescale, 5.1f*flarescale);
				const float flare10Bpow = 1.5f;
				const float flare10Bfill = 1.0f;
				const float flare10Boffset = -0.77f;
				vec2 flare10Bpos = vec2(  ((1.0 - lPos.x)*(flare10Boffset + 1.0) - (flare10Boffset*0.5))  *aspectRatio*flare10Bscale.x,  ((1.0 - lPos.y)*(flare10Boffset + 1.0) - (flare10Boffset*0.5))  *flare10Bscale.y);
			
			
				float flare10B = distance(flare10Bpos, vec2(texcoord.s*aspectRatio*flare10Bscale.x, texcoord.t*flare10Bscale.y));
					flare10B = 0.5 - flare10B;
					flare10B = clamp(flare10B*flare10Bfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
					flare10B = sin(flare10B*1.57075);
					flare10B *= sunmask;
					flare10B = pow(flare10B, 0.9f);
				  
					flare10B *= flare10Bpow;
				  
				flare10 = clamp(flare10 - flare10B, 0.0, 10.0);
				  
				  
				  	color.r += flare10*0.8f*flaremultR;
					color.g += flare10*0.2f*flaremultG;
					color.b += flare10*0.0f*flaremultB;				
					
					
			//mid blue sweep
			const vec2 flare10Cscale = vec2(6.0f*flarescale, 6.0f*flarescale);
			const float flare10Cpow = 1.9f;
			const float flare10Cfill = 1.1f;
			const float flare10Coffset = -0.6f;
			vec2 flare10Cpos = vec2(  ((1.0 - lPos.x)*(flare10Coffset + 1.0) - (flare10Coffset*0.5))  *aspectRatio*flare10Cscale.x,  ((1.0 - lPos.y)*(flare10Coffset + 1.0) - (flare10Coffset*0.5))  *flare10Cscale.y);
			
			
			float flare10C = distance(flare10Cpos, vec2(texcoord.s*aspectRatio*flare10Cscale.x, texcoord.t*flare10Cscale.y));
				  flare10C = 0.5 - flare10C;
				  flare10C = clamp(flare10C*flare10Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare10C = sin(flare10C*1.57075);
				  flare10C *= sunmask;
				  flare10C = pow(flare10C, 1.1f);
				  
				  flare10C *= flare10Cpow;			
				  
				  
				  //subtract
				const vec2 flare10Dscale = vec2(5.1f*flarescale, 5.1f*flarescale);
				const float flare10Dpow = 1.5f;
				const float flare10Dfill = 1.0f;
				const float flare10Doffset = -0.67f;
				vec2 flare10Dpos = vec2(  ((1.0 - lPos.x)*(flare10Doffset + 1.0) - (flare10Doffset*0.5))  *aspectRatio*flare10Dscale.x,  ((1.0 - lPos.y)*(flare10Doffset + 1.0) - (flare10Doffset*0.5))  *flare10Dscale.y);
			
			
				float flare10D = distance(flare10Dpos, vec2(texcoord.s*aspectRatio*flare10Dscale.x, texcoord.t*flare10Dscale.y));
					flare10D = 0.5 - flare10D;
					flare10D = clamp(flare10D*flare10Dfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
					flare10D = sin(flare10D*1.57075);
					flare10D *= sunmask;
					flare10D = pow(flare10D, 0.9f);
				  
					flare10D *= flare10Dpow;
				  
				flare10C = clamp(flare10C - flare10D, 0.0, 10.0);
				  
				  
				  	color.r += flare10C*0.0f*flaremultR;
					color.g += flare10C*0.2f*flaremultG;
					color.b += flare10C*0.9f*flaremultB;	
		//////////////////////////////////////////////////////////
		
		
		
		
		
		//Pointy fuzzy glow dots////////////////////////////////////////////////
			//RedGlow1

			const vec2 flare11scale = vec2(1.5f*flarescale, 1.5f*flarescale);
			const float flare11pow = 1.1f;
			const float flare11fill = 2.0f;
			const float flare11offset = -0.523f;
			vec2 flare11pos = vec2(  ((1.0 - lPos.x)*(flare11offset + 1.0) - (flare11offset*0.5))  *aspectRatio*flare11scale.x,  ((1.0 - lPos.y)*(flare11offset + 1.0) - (flare11offset*0.5))  *flare11scale.y);
			
			
			float flare11 = distance(flare11pos, vec2(texcoord.s*aspectRatio*flare11scale.x, texcoord.t*flare11scale.y));
				  flare11 = 0.5 - flare11;
				  flare11 = clamp(flare11*flare11fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare11 = pow(flare11, 2.9f);
				  flare11 *= sunmask;

				  
				  flare11 *= flare11pow;
				  
				  	color.r += flare11*1.0f*flaremultR;
					color.g += flare11*0.2f*flaremultG;
					color.b += flare11*0.0f*flaremultB;		
					
					
			//PurpleGlow2

			const vec2 flare12scale = vec2(2.5f*flarescale, 2.5f*flarescale);
			const float flare12pow = 0.5f;
			const float flare12fill = 2.0f;
			const float flare12offset = -0.323f;
			vec2 flare12pos = vec2(  ((1.0 - lPos.x)*(flare12offset + 1.0) - (flare12offset*0.5))  *aspectRatio*flare12scale.x,  ((1.0 - lPos.y)*(flare12offset + 1.0) - (flare12offset*0.5))  *flare12scale.y);
			
			
			float flare12 = distance(flare12pos, vec2(texcoord.s*aspectRatio*flare12scale.x, texcoord.t*flare12scale.y));
				  flare12 = 0.5 - flare12;
				  flare12 = clamp(flare12*flare12fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare12 = pow(flare12, 2.9f);
				  flare12 *= sunmask;

				  
				  flare12 *= flare12pow;
				  
				  	color.r += flare12*0.8f*flaremultR;
					color.g += flare12*0.2f*flaremultG;
					color.b += flare12*1.0f*flaremultB;		
					
					
					
			//BlueGlow3

			const vec2 flare13scale = vec2(1.0f*flarescale, 1.0f*flarescale);
			const float flare13pow = 0.5f;
			const float flare13fill = 2.0f;
			const float flare13offset = +0.138f;
			vec2 flare13pos = vec2(  ((1.0 - lPos.x)*(flare13offset + 1.0) - (flare13offset*0.5))  *aspectRatio*flare13scale.x,  ((1.0 - lPos.y)*(flare13offset + 1.0) - (flare13offset*0.5))  *flare13scale.y);
			
			
			float flare13 = distance(flare13pos, vec2(texcoord.s*aspectRatio*flare13scale.x, texcoord.t*flare13scale.y));
				  flare13 = 0.5 - flare13;
				  flare13 = clamp(flare13*flare13fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare13 = pow(flare13, 2.9f);
				  flare13 *= sunmask;

				  
				  flare13 *= flare13pow;
				  
				  	color.r += flare13*0.0f*flaremultR;
					color.g += flare13*0.2f*flaremultG;
					color.b += flare13*1.0f*flaremultB;		
					
			
			color.rgb = clamp(color.rgb, 0.0, 10.0);

#endif

/*
//AO FILTER

			float aosample = 0.0f;
			float weightsample = 0.0f;
			const float aospread = 11.2f * (1.0f - depth);
			const float aosharpness = 1.0f;
			
			
			
			
			float aosampleweights = 0.0f;
			
			float aosx = -0.030*aospread;
			float aosy = -0.030*aospread*aspectRatio;
			
			aosample += 1.0f - (1.0f - texture2D(gaux1, texcoord.st).a);
			
			for (int i = 0; i < 5; i++) {
			
					for (int i = 0; i < 5; i++) {
						weightsample = 1.0f - clamp((distance(getDepth(texcoord.st), getDepth(texcoord.st + vec2(aosx, aosy))))*5.0f, 0.0f, 1.0f);
						aosample += 1.0f - (1.0f - texture2D(gaux1, texcoord.st + vec2(aosx, aosy)).a * weightsample);
				
						aosy += 0.01*aospread;
						aosampleweights += weightsample;
					}

				aosx += 0.01*aospread;
				aosy = 0.01*aospread*aspectRatio;
				
				}
				
				aosample /= aosampleweights + 1.0f;
				
				color.rgb *= aosample;
*/

//color.rgb *= texture2D(gaux1, texcoord.st).a;


#ifdef CEL_SHADING
	color.rgb *= (getCellShaderFactor(texcoord.st));
#endif



#ifdef HDR


#endif;

color = color * BRIGHTMULT;

#ifdef CROSSPROCESS
	//pre-gain
	color = color * (BRIGHTMULT + 0.0f) + 0.03f;
	
	//compensate for low-light artifacts
	color = color+0.029f;
 
	//calculate double curve
	float dbr = -color.r + 1.4f;
	float dbg = -color.g + 1.4f;
	float dbb = -color.b + 1.4f;
	
	//fade between simple gamma up curve and double curve
	float pr = mix(dbr, 0.55f, 0.7f);
	float pg = mix(dbg, 0.55f, 0.7f);
	float pb = mix(dbb, 0.55f, 0.7f);
	
	color.r = pow((color.r * 0.95f - 0.005f), pr);
	color.g = pow((color.g * 0.95f - 0.002f), pg);
	color.b = pow((color.b * 0.91f + 0.000f), pb);
#endif

	//Color boosting
	color.r = (color.r)*(COLOR_BOOST + 1.0f) + (color.g + color.b)*(-COLOR_BOOST);
	color.g = (color.g)*(COLOR_BOOST + 1.0f) + (color.r + color.b)*(-COLOR_BOOST);
	color.b = (color.b)*(COLOR_BOOST + 1.0f) + (color.r + color.g)*(-COLOR_BOOST);
	
	//color.r = mix(((color.r)*(COLOR_BOOST + 1.0) + (hld.g + hld.b)*(-COLOR_BOOST)), hld.r, (max(((1-rgb)*2 - 1), 0.0)));
	//color.g = mix(((color.g)*(COLOR_BOOST + 1.0) + (hld.r + hld.b)*(-COLOR_BOOST)), hld.g, (max(((1-rgb)*2 - 1), 0.0)));
	//color.b = mix(((color.b)*(COLOR_BOOST + 1.0) + (hld.r + hld.g)*(-COLOR_BOOST)), hld.b, (max(((1-rgb)*2 - 1), 0.0)));

#ifdef HIGHDESATURATE


	//average
	float rgb = max(color.r, max(color.g, color.b))/2 + min(color.r, min(color.g, color.b))/2;

	//adjust black and white image to be brighter
	float bw = pow(rgb, 0.7f);

	//mix between per-channel analysis and average analysis
	float rgbr = mix(rgb, color.r, 0.7f);
	float rgbg = mix(rgb, color.g, 0.7f);
	float rgbb = mix(rgb, color.b, 0.7f);

	//calculate crossfade based on lum
	float mixfactorr = max(0.0f, (rgbr*4.0f - 3.0f));
	float mixfactorg = max(0.0f, (rgbg*4.0f - 3.0f));
	float mixfactorb = max(0.0f, (rgbb*4.0f - 3.0f));

	//crossfade between saturated and desaturated image
	float mixr = mix(color.r, bw, mixfactorr);
	float mixg = mix(color.g, bw, mixfactorg);
	float mixb = mix(color.b, bw, mixfactorb);

	//adjust level of desaturation
	color.r = clamp((mix(mixr, color.r, 1.0)), 0.0f, 10.0f);
	color.g = clamp((mix(mixg, color.g, 1.0)), 0.0f, 10.0f);
	color.b = clamp((mix(mixb, color.b, 1.0)), 0.0f, 10.0f);
	
	//desaturate blue channel
	//color.b = color.b*0.8f + ((color.r + color.g)/2.0f)*0.2f;
	

	//hold color values for color boost
	//vec4 hld = color;

	
	

	

	
	//color = color * BRIGHTMULT;


	
#endif

	//undo artifact compensation
	color = max(((color*1.10f) - 0.06f), 0.0f);

	color.r = pow(color.r, GAMMA);
	color.g = pow(color.g, GAMMA);
	color.b = pow(color.b, GAMMA);

	

	
//color *= 1.1f;

#ifdef VINTAGE

	color.r = clamp(color.r, 0.04, 1.0);

	color.b = clamp(color.b, 0.06, 0.89);
	
	//color.r = pow(color.r, GAMMA);
	//color.g = pow(color.g, GAMMA);
	//color.b = pow(color.b, GAMMA);
#endif



#ifdef TONEMAP

	/*
	float adaptation = 0.0f;
	float Yhdrsample = 0.20f;
	int HDRsamples = 0;
	const float HDRnoiseamp = 0.00f;
	noiseX4 = 0.0f;
	noiseY4 = 0.0f;
	
		for (int i = 0; i < 5; ++i) {
	
			adaptation += texture2D(gaux1, vec2(0.30f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.35f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.40f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.45f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.50f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.55f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.60f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.65f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.70f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.75f, Yhdrsample) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			
			
			//adaptation += texture2D(gaux1, vec2(0.305f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.355f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.405f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.455f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.505f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.555f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.605f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.655f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			//adaptation += texture2D(gaux1, vec2(0.705f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			adaptation += texture2D(gaux1, vec2(0.755f, Yhdrsample+0.005) + vec2(noiseX4*HDRnoiseamp, noiseY4*HDRnoiseamp)).a;
			
			
			Yhdrsample += 0.05f;
			HDRsamples += 10;

		}
		
		adaptation /= HDRsamples;
		float adaptation_grayscale = adaptation;
				//adaptation_grayscale = pow(adaptation_grayscale, 1.1f);
				adaptation_grayscale = clamp(adaptation_grayscale, 0.0f, 1.0f);
				adaptation_grayscale *= 0.4f;
				adaptation_grayscale += 0.12f;
		
		
		
	const float interiorwarmth = 0.3f;
	const float toneboost = 1.6f;

		
	//Exposure
	color.r = color.r * toneboost / (adaptation_grayscale);
	color.g = color.g * toneboost / (mix(adaptation_grayscale, adaptation_grayscale * 0.8f + 0.1f, interiorwarmth));
	color.b = color.b * toneboost / (mix(adaptation_grayscale, adaptation_grayscale * 0.5f + 0.25f, interiorwarmth));
	color.b *= 1.1f;
*/



	color.rgb *= 1.1f;
	const float TonemapOversat = 7.0f;
	const float TonemapCurve = 0.28;

	//Tonemap
	color.rgb = (color.rgb * (1.0 + color.rgb/TonemapOversat))/(color.rgb + TonemapCurve);
	
	color = color*(1.0f + DARKMULT) - DARKMULT;
	
	color = clamp(color, 0.0f, 1.0f);

#endif



	gl_FragColor = color;
	
// End of Main. -----------------
}
