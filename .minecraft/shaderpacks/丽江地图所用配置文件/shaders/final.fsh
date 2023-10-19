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




//#define BOKEH_DOF							//Cannot be applied to water
#define HQ_DOF								//Enable for higher quality DOF
#define TILT_SHIFT							//Tilt shift effect. Not meant for gameplay. Google "tilt shift" for more info.
#define TILT_SHIFT_SCALE 0.5				//Size of aperture. Higher values gives illusion of smaller world
#define GODRAYS
#define GODRAYS_EXPOSURE 0.10
#define GODRAYS_SAMPLES 6
#define GODRAYS_DECAY 0.99
#define GODRAYS_DENSITY 0.30
//#define LENS								//ATI Cards only
#define LENS_POWER 0.56
//#define GLARE
#define GLARE_AMOUNT 0.25
#define GLARE_RANGE 2.0
//#define GLARE2							//second pass of glare shader. More realistic light scattering.
//#define CEL_SHADING
//#define CEL_SHADING_THRESHOLD 0.4
//#define CEL_SHADING_THICKNESS 0.004




//#define VINTAGE
#define VIGNETTE
#define VIGNETTE_STRENGTH 1.30
#define LOWLIGHT_EYE
#define TONEMAP	
#define TONEMAP_FILMIC						
#define TONEMAP_COLOR_FILTER
#define BRIGHTMULT 1.00               	// 1.0 = default brightness. Higher values mean brighter. 0 would be black.
#define DARKMULT 0.00					// 0.0 = normal image. Higher values will darken dark colors.

#ifdef TONEMAP_FILMIC
#define COLOR_BOOST	0.09
#else
#define COLOR_BOOST	0.13				// 0.0 = normal saturation. Higher values mean more saturated image.
#endif

//#define MOTIONBLUR					// Cannot be applied to water
#define MOTIONBLUR_AMOUNT 1.5
#define GAMMA 1.00f						//1.0 is default brightness. lower values will brighten image, higher values will darken image	

#define WATER_SHADER
#define GLOSSY_REFLECTIONS

//#define SCREEN_SPACE_RADIOSITY
#define RADIOSITY_AMOUNT 1.1

#define ICE_SHADER





// DOF Constants - DO NOT CHANGE
// HYPERFOCAL = (Focal Distance ^ 2)/(Circle of Confusion * F Stop) + Focal Distance
#ifdef USE_DOF
const float HYPERFOCAL = 3.132;
const float PICONSTANT = 3.14159;
#endif





//uniform sampler2D texture;
uniform sampler2D gdepth;
uniform sampler2D gcolor;
uniform sampler2D gnormal;
uniform sampler2D composite;
uniform sampler2D gaux1; // red is our motion blur mask. If red == 1, don't blur. green is water mask
uniform sampler2D gaux2; // red is godrays
uniform sampler2D gaux3; 

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferPreviousProjection;

uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferPreviousModelView;

uniform vec3 cameraPosition;
uniform vec3 previousCameraPosition;

uniform vec3 sunPosition;

uniform int worldTime;
uniform float aspectRatio;
uniform float near;
uniform float far;
uniform float viewWidth;
uniform float viewHeight;
uniform float rainStrength;
uniform float wetness;

varying vec4 texcoord;



//Land/sky mask
float land = texture2D(gaux1, texcoord.st).b;

//Raining
float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;
float wetx  = clamp(wetness, 0.0f, 1.0f);

//Water mask
//float iswater = texture2D(gaux1, texcoord.st).g;
float isice   = texture2D(gaux3, texcoord.st).b;

vec3 specularity = texture2D(gaux3, texcoord.st).rgb;


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

vec3 albedo = texture2D(gcolor, texcoord.st).rgb;

	
#ifdef BOKEH_DOF

	//compare dof depth function
	float dofWeight(vec2 blur, vec2 coord) {
		float dthresh = 500.0;
		float dthresh2 = 1.0f;
		return (1.0f - (clamp((texture2D(gdepth, texcoord.st).x - texture2D(gdepth, texcoord.st + coord).x) * dthresh, 0.0f, 1.0f)) * (1.0f - clamp(abs(blur.x) * dthresh2, 0.0f, 1.0f)));	
		//return 1.0f;	
	}

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
			
			for (int i = 0; i < 2; ++i) {
			
			if (texcoord.s > 1.0f || texcoord.s < 0.0f || texcoord.t > 1.0f || texcoord.t < 0.0f) {
				break;
			}
			
				
				float sample = 0.0f;

					sample = 1.0f - texture2D(gaux2, tx + vec2(noise*delta.x, noise*delta.y)).r;
					sample += 1.0f - texture2D(gaux2, tx + vec2(noise2*delta.x, noise2*delta.y)).r;
					sample += 1.0f - texture2D(gaux2, tx + vec2(noise3*delta.x, noise3*delta.y)).r;
					sample += 1.0f - texture2D(gaux2, tx + vec2(noise4*delta.x, noise4*delta.y)).r;
					sample += 1.0f - texture2D(gaux2, tx + vec2(noise5*delta.x, noise5*delta.y)).r;
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
	
	int iswater;
	
	if(texture2D(gaux1, texcoord.st).g == 1.0f){
		iswater = 1;
	} else {
		iswater = 0;
	}
	
	
//Common variables

	float depth = eDepth(texcoord.xy);
	vec2 Texcoord2 = texcoord.st;
	float linDepth = getDepth(texcoord.st);
	vec3 normal = texture2D(gnormal, texcoord.st).rgb;
	vec3 normalBiased = normal * 2.0f - 1.0f;

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
					
						float width7 = 6.0;
						float height7 = 6.0;
						float noiseX7 = ((fract(1.0-Texcoord2.s*(width7/2.0))*0.25)+(fract(Texcoord2.t*(height7/2.0))*0.75))*2.0-1.0;
						float noiseY7 = ((fract(1.0-Texcoord2.s*(width7/2.0))*0.75)+(fract(Texcoord2.t*(height7/2.0))*0.25))*2.0-1.0;

						
							noiseX7 = clamp(fract(sin(dot(Texcoord2 ,vec2(12.9898,44.633))) * 51178.5453),0.0,1.0)*2.0-1.0;
							noiseY7 = clamp(fract(sin(dot(Texcoord2 ,vec2(43.9898,61.233)*2.0)) * 9958.5453),0.0,1.0)*2.0-1.0;
						
						noiseX7 *= (0.10f*noiseamp);
						noiseY7 *= (0.10f*noiseamp);
						
						float width8 = 7.0;
						float height8 = 7.0;
						float noiseX8 = ((fract(1.0-Texcoord2.s*(width8/2.0))*0.25)+(fract(Texcoord2.t*(height8/2.0))*0.75))*2.0-1.0;
						float noiseY8 = ((fract(1.0-Texcoord2.s*(width8/2.0))*0.75)+(fract(Texcoord2.t*(height8/2.0))*0.25))*2.0-1.0;

						
							noiseX8 = clamp(fract(sin(dot(Texcoord2 ,vec2(14.9898,47.633))) * 51468.5453),0.0,1.0)*2.0-1.0;
							noiseY8 = clamp(fract(sin(dot(Texcoord2 ,vec2(13.9898,81.233)*2.0)) * 6388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX8 *= (0.10f*noiseamp);
						noiseY8 *= (0.10f*noiseamp);
						
						float width9 = 8.0;
						float height9 = 8.0;
						float noiseX9 = ((fract(1.0-Texcoord2.s*(width9/2.0))*0.25)+(fract(Texcoord2.t*(height9/2.0))*0.75))*2.0-1.0;
						float noiseY9 = ((fract(1.0-Texcoord2.s*(width9/2.0))*0.75)+(fract(Texcoord2.t*(height9/2.0))*0.25))*2.0-1.0;

						
							noiseX9 = clamp(fract(sin(dot(Texcoord2 ,vec2(24.9898,59.633))) * 55468.5453),0.0,1.0)*2.0-1.0;
							noiseY9 = clamp(fract(sin(dot(Texcoord2 ,vec2(23.9898,95.233)*2.0)) * 16388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX9 *= (0.10f*noiseamp);
						noiseY9 *= (0.10f*noiseamp);
						
						float width10 = 9.0;
						float height10 = 9.0;
						float noiseX10 = ((fract(1.0-Texcoord2.s*(width10/2.0))*0.25)+(fract(Texcoord2.t*(height10/2.0))*0.75))*2.0-1.0;
						float noiseY10 = ((fract(1.0-Texcoord2.s*(width10/2.0))*0.75)+(fract(Texcoord2.t*(height10/2.0))*0.25))*2.0-1.0;

						
							noiseX10 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898,59.633))) * 57468.5453),0.0,1.0)*2.0-1.0;
							noiseY10 = clamp(fract(sin(dot(Texcoord2 ,vec2(25.9898,95.233)*2.0)) * 18388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX10 *= (0.10f*noiseamp);
						noiseY10 *= (0.10f*noiseamp);
					
					/*
						float width11 = 10.0;
						float height11 = 10.0;
						float noiseX11 = ((fract(1.0-Texcoord2.s*(width11/2.0))*0.25)+(fract(Texcoord2.t*(height11/2.0))*0.75))*2.0-1.0;
						float noiseY11 = ((fract(1.0-Texcoord2.s*(width11/2.0))*0.75)+(fract(Texcoord2.t*(height11/2.0))*0.25))*2.0-1.0;

						
							noiseX11 = clamp(fract(sin(dot(Texcoord2 ,vec2(28.9898,61.633))) * 59468.5453),0.0,1.0)*2.0-1.0;
							noiseY11 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898,97.233)*2.0)) * 21388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX11 *= 0.002;
						noiseY11 *= 0.002;
						
						float width12 = 11.0;
						float height12 = 11.0;
						float noiseX12 = ((fract(1.0-Texcoord2.s*(width12/2.0))*0.25)+(fract(Texcoord2.t*(height12/2.0))*0.75))*2.0-1.0;
						float noiseY12 = ((fract(1.0-Texcoord2.s*(width12/2.0))*0.75)+(fract(Texcoord2.t*(height12/2.0))*0.25))*2.0-1.0;

						
							noiseX12 = clamp(fract(sin(dot(Texcoord2 ,vec2(30.9898,64.633))) * 61468.5453),0.0,1.0)*2.0-1.0;
							noiseY12 = clamp(fract(sin(dot(Texcoord2 ,vec2(34.9898,99.233)*2.0)) * 23388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX12 *= 0.002;
						noiseY12 *= 0.002;		

						float width13 = 12.0;
						float height13 = 12.0;
						float noiseX13 = ((fract(1.0-Texcoord2.s*(width13/2.0))*0.25)+(fract(Texcoord2.t*(height13/2.0))*0.75))*2.0-1.0;
						float noiseY13 = ((fract(1.0-Texcoord2.s*(width13/2.0))*0.75)+(fract(Texcoord2.t*(height13/2.0))*0.25))*2.0-1.0;

						
							noiseX13 = clamp(fract(sin(dot(Texcoord2 ,vec2(32.9898,66.633))) * 63468.5453),0.0,1.0)*2.0-1.0;
							noiseY13 = clamp(fract(sin(dot(Texcoord2 ,vec2(36.9898,101.233)*2.0)) * 25388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX13 *= 0.002;
						noiseY13 *= 0.002;		

						float width14 = 13.0;
						float height14 = 13.0;
						float noiseX14 = ((fract(1.0-Texcoord2.s*(width14/2.0))*0.25)+(fract(Texcoord2.t*(height14/2.0))*0.75))*2.0-1.0;
						float noiseY14 = ((fract(1.0-Texcoord2.s*(width14/2.0))*0.75)+(fract(Texcoord2.t*(height14/2.0))*0.25))*2.0-1.0;

						
							noiseX14 = clamp(fract(sin(dot(Texcoord2 ,vec2(34.9898,68.633))) * 65468.5453),0.0,1.0)*2.0-1.0;
							noiseY14 = clamp(fract(sin(dot(Texcoord2 ,vec2(38.9898,103.233)*2.0)) * 27388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX14 *= 0.002;
						noiseY14 *= 0.002;	

						float width15 = 14.0;
						float height15 = 14.0;
						float noiseX15 = ((fract(1.0-Texcoord2.s*(width15/2.0))*0.25)+(fract(Texcoord2.t*(height15/2.0))*0.75))*2.0-1.0;
						float noiseY15 = ((fract(1.0-Texcoord2.s*(width15/2.0))*0.75)+(fract(Texcoord2.t*(height15/2.0))*0.25))*2.0-1.0;

						
							noiseX15 = clamp(fract(sin(dot(Texcoord2 ,vec2(36.9898,70.633))) * 67468.5453),0.0,1.0)*2.0-1.0;
							noiseY15 = clamp(fract(sin(dot(Texcoord2 ,vec2(40.9898,105.233)*2.0)) * 29388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX15 *= 0.002;
						noiseY15 *= 0.002;	

						float width16 = 15.0;
						float height16 = 15.0;
						float noiseX16 = ((fract(1.0-Texcoord2.s*(width16/2.0))*0.25)+(fract(Texcoord2.t*(height16/2.0))*0.75))*2.0-1.0;
						float noiseY16 = ((fract(1.0-Texcoord2.s*(width16/2.0))*0.75)+(fract(Texcoord2.t*(height16/2.0))*0.25))*2.0-1.0;

						
							noiseX16 = clamp(fract(sin(dot(Texcoord2 ,vec2(38.9898,72.633))) * 69468.5453),0.0,1.0)*2.0-1.0;
							noiseY16 = clamp(fract(sin(dot(Texcoord2 ,vec2(42.9898,107.233)*2.0)) * 31388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX16 *= 0.002;
						noiseY16 *= 0.002;		

						float width17 = 16.0;
						float height17 = 16.0;
						float noiseX17 = ((fract(1.0-Texcoord2.s*(width17/2.0))*0.25)+(fract(Texcoord2.t*(height17/2.0))*0.75))*2.0-1.0;
						float noiseY17 = ((fract(1.0-Texcoord2.s*(width17/2.0))*0.75)+(fract(Texcoord2.t*(height17/2.0))*0.25))*2.0-1.0;

						
							noiseX17 = clamp(fract(sin(dot(Texcoord2 ,vec2(40.9898,74.633))) * 70468.5453),0.0,1.0)*2.0-1.0;
							noiseY17 = clamp(fract(sin(dot(Texcoord2 ,vec2(44.9898,109.233)*2.0)) * 33388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX17 *= 0.002;
						noiseY17 *= 0.002;

						float width18 = 17.0;
						float height18 = 17.0;
						float noiseX18 = ((fract(1.0-Texcoord2.s*(width18/2.0))*0.25)+(fract(Texcoord2.t*(height18/2.0))*0.75))*2.0-1.0;
						float noiseY18 = ((fract(1.0-Texcoord2.s*(width18/2.0))*0.75)+(fract(Texcoord2.t*(height18/2.0))*0.25))*2.0-1.0;

						
							noiseX18 = clamp(fract(sin(dot(Texcoord2 ,vec2(42.9898,76.633))) * 72468.5453),0.0,1.0)*2.0-1.0;
							noiseY18 = clamp(fract(sin(dot(Texcoord2 ,vec2(46.9898,111.233)*2.0)) * 35388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX18 *= 0.002;
						noiseY18 *= 0.002;	

						float width19 = 18.0;
						float height19 = 18.0;
						float noiseX19 = ((fract(1.0-Texcoord2.s*(width19/2.0))*0.25)+(fract(Texcoord2.t*(height19/2.0))*0.75))*2.0-1.0;
						float noiseY19 = ((fract(1.0-Texcoord2.s*(width19/2.0))*0.75)+(fract(Texcoord2.t*(height19/2.0))*0.25))*2.0-1.0;

						
							noiseX19 = clamp(fract(sin(dot(Texcoord2 ,vec2(44.9898,78.633))) * 75468.5453),0.0,1.0)*2.0-1.0;
							noiseY19 = clamp(fract(sin(dot(Texcoord2 ,vec2(48.9898,115.233)*2.0)) * 38388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX19 *= 0.002;
						noiseY19 *= 0.002;		

						float width20 = 19.0;
						float height20 = 19.0;
						float noiseX20 = ((fract(1.0-Texcoord2.s*(width20/2.0))*0.25)+(fract(Texcoord2.t*(height20/2.0))*0.75))*2.0-1.0;
						float noiseY20 = ((fract(1.0-Texcoord2.s*(width20/2.0))*0.75)+(fract(Texcoord2.t*(height20/2.0))*0.25))*2.0-1.0;

						
							noiseX20 = clamp(fract(sin(dot(Texcoord2 ,vec2(46.9898,81.633))) * 77468.5453),0.0,1.0)*2.0-1.0;
							noiseY20 = clamp(fract(sin(dot(Texcoord2 ,vec2(51.9898,118.233)*2.0)) * 41188.5453),0.0,1.0)*2.0-1.0;
						
						noiseX20 *= 0.002;
						noiseY20 *= 0.002;		

						float width21 = 20.0;
						float height21 = 20.0;
						float noiseX21 = ((fract(1.0-Texcoord2.s*(width21/2.0))*0.25)+(fract(Texcoord2.t*(height21/2.0))*0.75))*2.0-1.0;
						float noiseY21 = ((fract(1.0-Texcoord2.s*(width21/2.0))*0.75)+(fract(Texcoord2.t*(height21/2.0))*0.25))*2.0-1.0;

						
							noiseX21 = clamp(fract(sin(dot(Texcoord2 ,vec2(48.9898,83.633))) * 79468.5453),0.0,1.0)*2.0-1.0;
							noiseY21 = clamp(fract(sin(dot(Texcoord2 ,vec2(53.9898,120.233)*2.0)) * 43188.5453),0.0,1.0)*2.0-1.0;
						
						noiseX21 *= 0.002;
						noiseY21 *= 0.002;	

						float width22 = 21.0;
						float height22 = 21.0;
						float noiseX22 = ((fract(1.0-Texcoord2.s*(width22/2.0))*0.25)+(fract(Texcoord2.t*(height22/2.0))*0.75))*2.0-1.0;
						float noiseY22 = ((fract(1.0-Texcoord2.s*(width22/2.0))*0.75)+(fract(Texcoord2.t*(height22/2.0))*0.25))*2.0-1.0;

						
							noiseX22 = clamp(fract(sin(dot(Texcoord2 ,vec2(51.9898,83.633))) * 81468.5453),0.0,1.0)*2.0-1.0;
							noiseY22 = clamp(fract(sin(dot(Texcoord2 ,vec2(56.9898,120.233)*2.0)) * 48188.5453),0.0,1.0)*2.0-1.0;
						
						noiseX22 *= 0.002;
						noiseY22 *= 0.002;
*/
#ifdef BOKEH_DOF
	
	if (depth > 0.9999f) {
		depth = 1.0f;
	}
	

	float cursorDepth = eDepth(vec2(0.5f, 0.5f));
	
	if (cursorDepth > 0.9999f) {
		cursorDepth = 1.0f;
	}
	
float blurclamp = 0.014;  // max blur amount
float bias = 0.3;	//aperture - bigger values for shallower depth of field

#ifdef TILT_SHIFT

	bias *= 80.0f * TILT_SHIFT_SCALE;
	blurclamp *= 80.0f * TILT_SHIFT_SCALE;

#endif
	
	
	vec2 aspectcorrect = vec2(1.0, aspectRatio) * 1.5;
	
	float factor = (depth - cursorDepth);
	 
	vec2 dofblur = (vec2 (clamp( factor * bias, -blurclamp, blurclamp )))*0.6;

	
	#ifdef HQ_DOF
	


	//HQ
	vec3 col = vec3(0.0);
	float dweight;
	float dweightall;
	
	col += texture2D(composite, texcoord.st).rgb;
						  
						  
						dweight =   dofWeight(dofblur, (vec2( 0.0, 0.4)*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0, 0.4)*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( 0.15,0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur).rgb * dweight;
						
						dweight =   dofWeight(dofblur, (vec2( 0.29,0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( -0.37,0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(composite, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	

	color.rgb = col/(dweightall + 0.001);	
	
	
	
	#else
	
	//LQ
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
	
#endif


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
		
		//float velocityadd = texture2D(gaux1, texcoord.st).a;
		//velocity += velocityadd;
		
		
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
			
			coord = texcoord.st + (velocity);

			//color += texture2D(composite, coord - vec2(noiseX2*velocity.x, noiseX2*velocity.y));
			//color += texture2D(composite, coord - vec2(noiseY2*velocity.x, noiseY2*velocity.y));
			//color += texture2D(composite, coord - vec2(noiseX4*velocity.x, noiseX4*velocity.y));
			//color += texture2D(composite, coord - vec2(noiseY4*velocity.x, noiseY4*velocity.y));
			color = texture2D(composite, texcoord.st + velocity);
			samples += 1;
			
			offsetcount += 1;
			
			coord = texcoord.st;
		
		}
			color = (color/1.0)/samples;
		}
		
		

	
#endif


/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
#ifdef WATER_SHADER
	
	const float rspread = 0.30f;						//How long reflections are spread across the screen
	
	float rdepth = depth;


float pix_x = 1.0f / viewWidth;
float pix_y = 1.0f / viewHeight;

	rdepth = pow(rdepth, 1.0f);

const float wnormalclamp = 0.05f;
	
//Detect water surface normals

	//Compare change in depth texture over 1 pixel and return an angle
		float wnormal_x1 = texture2D(gdepth, texcoord.st + vec2(pix_x, 0.0f)).x - texture2D(gdepth, texcoord.st).x;
		float wnormal_x2 = texture2D(gdepth, texcoord.st).x - texture2D(gdepth, texcoord.st + vec2(-pix_x, 0.0f)).x;			
		float wnormal_x = 0.0f;
		
		if(abs(wnormal_x1) > abs(wnormal_x2)){
			wnormal_x = wnormal_x2;
		} else {
			wnormal_x = wnormal_x1;
		}
		wnormal_x /= 1.0f - rdepth;	

		wnormal_x = clamp(wnormal_x, -wnormalclamp, wnormalclamp);
		
		wnormal_x *= rspread*1.0f;
		

			  
			  
		float wnormal_y1 = texture2D(gdepth, texcoord.st + vec2(0.0f, pix_y)).x - texture2D(gdepth, texcoord.st).x;
		float wnormal_y2 = texture2D(gdepth, texcoord.st).x - texture2D(gdepth, texcoord.st + vec2(0.0f, -pix_y)).x;		
		float wnormal_y;
		
		if(abs(wnormal_y1) > abs(wnormal_y2)){
			wnormal_y = wnormal_y2;
		} else {
			wnormal_y = wnormal_y1;
		}	
		wnormal_y /= 1.0f - rdepth;			

		wnormal_y = clamp(wnormal_y, -wnormalclamp, wnormalclamp);
		
		wnormal_y *= rspread*1.0f*aspectRatio;
		
		  		  
		  
		 
		  
		  
		//if (down >= 1.0f) {
		//		down = 0.0f;
		// }
          
		  
//REFRACTION

	//Heightmap of small waves
	float waves = texture2D(gaux2, texcoord.st).g;
	float wavesraw = waves;
		  waves -= 0.5f;
		  waves *= 1.0 - depth;
		  waves *= 100.0f;

	//Detect angle of waves by comparing 1 pixel difference and resolving discontinuities
	float wavesdeltax1 = texture2D(gaux2, texcoord.st).g - texture2D(gaux2, texcoord.st + vec2(-pix_x, 0.0f)).g;
	float wavesdeltax2 = texture2D(gaux2, texcoord.st + vec2(pix_x, 0.0f)).g - texture2D(gaux2, texcoord.st).g;
	float wavesdeltax;
	
		if(abs(wavesdeltax1) > abs(wavesdeltax2)){
			wavesdeltax = wavesdeltax2;
		} else {
			wavesdeltax = wavesdeltax1;
		}
		
		wavesdeltax = clamp(wavesdeltax, -0.1f, 0.1f);
		
		wavesdeltax *= 1.0f - depth;
		wavesdeltax *= 30.0f;
		  
		  
	float wavesdeltay1 = texture2D(gaux2, texcoord.st).g - texture2D(gaux2, texcoord.st + vec2(0.0f, -pix_y)).g;
	float wavesdeltay2 = texture2D(gaux2, texcoord.st + vec2(0.0f, pix_y)).g - texture2D(gaux2, texcoord.st).g;
	float wavesdeltay = 0.0f;
	
		if(abs(wavesdeltay1) > abs(wavesdeltay2)){
			wavesdeltay = wavesdeltay2;
		} else {
			wavesdeltay = wavesdeltay1;
		}
		wavesdeltay *= 1.0f - depth;
		wavesdeltay *= 30.0f;
		
		wavesdeltay = clamp(wavesdeltay, -0.1f, 0.1f);
		  


	


float refractamount = 500.1154f*1.75f;
float refractamount2 = 0.0214f*0.00f;
float refractamount3 = 0.214f*0.25f;
float waberration = 0.10f;

	vec3 refracted = vec3(0.0f);
	vec3 refractedmask = vec3(0.0f);
	float bigWaveRefract = 1000.0f * (1.0f - depth);
	float bigWaveRefractScale = 1500.0f * (1.0f - depth);
	
	vec2 bigRefract = vec2(wnormal_x*bigWaveRefract, wnormal_y*bigWaveRefract);
	
	for (int i = 0; i < 1; ++i) {
			
				if(iswater != 1.0f) {
					break;
				}
	
			vec2 refractcoord_r = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y) * (waberration * 2.0f + 1.0f);
			vec2 refractcoord_g = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y) * (waberration + 1.0f);
			vec2 refractcoord_b = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y);
				
				refractcoord_r = refractcoord_r * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				refractcoord_g = refractcoord_g * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				refractcoord_b = refractcoord_b * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				
			/*
			refractcoord_r.s = clamp(refractcoord_r.s, 0.001f, 0.999f);
			refractcoord_r.t = clamp(refractcoord_r.t, 0.001f, 0.999f);	
			
			refractcoord_g.s = clamp(refractcoord_g.s, 0.001f, 0.999f);
			refractcoord_g.t = clamp(refractcoord_g.t, 0.001f, 0.999f);
			
			refractcoord_b.s = clamp(refractcoord_b.s, 0.001f, 0.999f);
			refractcoord_b.t = clamp(refractcoord_b.t, 0.001f, 0.999f);
			*/
			
			
			if (refractcoord_r.s > 1.0 || refractcoord_r.s < 0.0 || refractcoord_r.t > 1.0 || refractcoord_r.t < 0.0 ||
				refractcoord_g.s > 1.0 || refractcoord_g.s < 0.0 || refractcoord_g.t > 1.0 || refractcoord_g.t < 0.0 ||
				refractcoord_b.s > 1.0 || refractcoord_b.s < 0.0 || refractcoord_b.t > 1.0 || refractcoord_b.t < 0.0) {
					break;
				}
			
			if (refractcoord_r.st == vec2(0.0f)) {
				break;
			}			
			
			if (refractcoord_g.st == vec2(0.0f)) {
				break;
			}			
			
			if (refractcoord_b.st == vec2(0.0f)) {
				break;
			}
			
			refracted.r = texture2D(composite, refractcoord_r).r;
			refracted.g = texture2D(composite, refractcoord_g).g;
			refracted.b = texture2D(composite, refractcoord_b).b;
			
			
			refractedmask.r = texture2D(gaux1, refractcoord_r).g;
			refractedmask.g = texture2D(gaux1, refractcoord_g).g;
			refractedmask.b = texture2D(gaux1, refractcoord_b).g;
	
			}
			
	color.r = mix(color.r, refracted.r, refractedmask.r);
	color.g = mix(color.g, refracted.g, refractedmask.g);
	color.b = mix(color.b, refracted.b, refractedmask.b);


//REFLECTION


	//color.rgb = worldposition.g;
	
	vec3 reflection = vec3(0.0f);
	float rtransy = 0.01f * rspread;
	float rtransin = 0.05f;
	
	const float rstrong = 5.4f;
	const float reflectwaviness = 0.00395f;
	const float rcurve = 1.0f;
	
	//coordinates for translating reflection
	vec2 coordnormal = vec2(0.0f);
	vec2 coordin = texcoord.st;
	vec2 rcoord = vec2(0.0f);
	
	float dwaves = waves * 0.4f * reflectwaviness;
	float dwavesdeltax = wavesdeltax * 7.3f * reflectwaviness;
	float dwavesdeltay = wavesdeltay * 7.3f * reflectwaviness;
	float reflectmask = 0.0f;
	float reflectmaskhold = 0.0f;
	float rnoise = 0.0f;
	
	float depthcheck = 0.0f;
	float depthcheck2 = 0.0f;
	float depthpass = 0.0f;
	float prevdepth;
	float thisdepth;
	
	int samples = 1;
	
	
			float redge = distance(texcoord.s, 0.5f);
			  redge = max(redge, distance(texcoord.t, 0.5f));
			  redge *= 2.0f;
			  redge = clamp(redge * 4.0f - 3.0f, 0.0f, 1.0f);
			  redge = 1.0f;
	
			
			for (int i = 0; i < 8; ++i) {
			
				if(iswater != 1.0f) {
					samples += 1;
					break;
				}
				
				rcoord = coordnormal + vec2(dwavesdeltax*4.0f + wnormal_x, dwavesdeltay*4.0f + wnormal_y)*(samples * samples - 1)*redge;
				
				thisdepth = texture2D(gdepth, clamp(texcoord.st + rcoord, 0.001f, 0.999f)).x;
				
				depthcheck = (rdepth - thisdepth);
				depthcheck = 1.0f - depthcheck;
				depthcheck = clamp(depthcheck * 140.0 - 139.0f, 0.0f, 1.0f);
				depthcheck2 = clamp(depthcheck * 70.0 - 69.0f, 0.0f, 1.0f);
				
				reflectmask   = ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseX2*wnormal_x*rnoise*(samples - 1), noiseX2*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).g) * ((9 - samples)/9.0f))/1.0f;
				//reflectmask  += ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseY2*wnormal_x*rnoise*(samples - 1), noiseY2*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).g) * ((9 - samples)/9.0f))/4.0f;
				//reflectmask  += ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseX4*wnormal_x*rnoise*(samples - 1), noiseX4*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).g) * ((9 - samples)/9.0f))/4.0f;
				//reflectmask  += ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseY4*wnormal_x*rnoise*(samples - 1), noiseY4*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).g) * ((9 - samples)/9.0f))/4.0f;
																																																																																																							
				reflection  += 	((texture2D(composite, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseX2*wnormal_x*rnoise*(samples - 1), noiseX2*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).rgb * reflectmask) * ((9 - samples)/9.0f))/1.0f;
				//reflection  += 	((texture2D(composite, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseY2*wnormal_x*rnoise*(samples - 1), noiseY2*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).rgb * reflectmask) * ((9 - samples)/9.0f))/4.0f;
				//reflection  += 	((texture2D(composite, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseX4*wnormal_x*rnoise*(samples - 1), noiseX4*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).rgb * reflectmask) * ((9 - samples)/9.0f))/4.0f;
				//reflection  += 	((texture2D(composite, clamp(texcoord.st + (rcoord*depthcheck + vec2(noiseY4*wnormal_x*rnoise*(samples - 1), noiseY4*wnormal_y*rnoise*(samples - 1))*redge*depthcheck), 0.001f, 0.999f)).rgb * reflectmask) * ((9 - samples)/9.0f))/4.0f;
				
				reflectmaskhold += reflectmask;

				samples += 1;
				

			
			}
			
			reflection /= samples - 1;
			reflectmaskhold /= samples - 1;
			
			reflectmaskhold = pow(reflectmaskhold, 1.0f)*2.5f;
			
			float wfresnel = pow(distance(vec2(wnormal_x, wnormal_y) + vec2(dwavesdeltax, dwavesdeltay), vec2(0.0f)), 0.7f) * 20.0f;
			
						
			//Darken objects behind water
			color.rgb = mix(color.rgb, vec3(color.r * (1.1f - wfresnel), color.g * (1.1f - wfresnel * 0.9f), color.b * (1.1f - wfresnel * 0.8f)) * (1.0 - reflectmaskhold), iswater);
			//color.rgb = mix(color.rgb, vec3(color.r * 0.05f, color.g * 0.05f, color.b * 0.05f) * (1.0 - reflectmaskhold), iswater);
			
			//Add reflections to water only >:3
			reflection *= iswater;
			
			color.rgb = color.rgb + (reflection * rstrong);
			
		//	float fakecaustic = abs(dwavesdeltax) + abs(dwavesdeltay);
		//	      fakecaustic *= 1000.0f;
		//		  fakecaustic = pow(fakecaustic, 10.0f);
		//		  fakecaustic *= iswater;
		//		  fakecaustic /= distance(vec2(-0.1), vec2(wnormal_x, wnormal_y));
				  
			//color.rgb *= 1.0f + max(fakecaustic, 0.0);
			
			
			

//faker reflections

		//color.rgb = mix(color.rgb, vec3(color.r * 0.25f, color.g * 0.35f, color.b * 0.35f), iswater);


			
#endif



///////////////ICE SHADER///////////////////
///////////////ICE SHADER///////////////////
///////////////ICE SHADER///////////////////
///////////////ICE SHADER///////////////////


#ifdef ICE_SHADER

//Detect ice surface normals
	float normalCurveI = 1.0f;
	
	
	float inormalX = texture2D(gdepth, texcoord.st + vec2(pix_x, 0.0f)).x - texture2D(gdepth, texcoord.st).x;
	      //inormalX = min(inormalX, texture2D(gdepth, texcoord.st + vec2(-pix_x, 0.0f).x - texture2D(gdepth, texcoord.st).x));
		
			if (inormalX < 0.0f) {
				inormalX = -inormalX;
				inormalX = pow(inormalX, normalCurveI);
				inormalX = -inormalX;
			}	else {
				inormalX = pow(inormalX, normalCurveI);
			}
			
		  inormalX /= 1.0f - depth;
			
		  inormalX *= 1.0f;
		  
		  
	float inormalY = texture2D(gdepth, texcoord.st + vec2(0.0f, pix_y)).x - texture2D(gdepth, texcoord.st).x;
			
			if (inormalY < 0.0f) {
				inormalY = -inormalY;
				inormalY = pow(inormalY, normalCurveI-0.00);
				inormalY = -inormalY;
			}	else {
				inormalY = pow(inormalY, normalCurveI-0.00);
			}
			
          inormalY /= 1.0f - depth;
			
		  inormalY *= 1.0f*aspectRatio;
		  
	
	
//refract//
/*
	vec3 IceRefract(float normalx, float normaly){
	
		if (isice != 1.0f) {
			break;
		}
	
		float refractStrength = 2000.0f * (1.0f - depth);
		float abb = 0.1f;
		
		vec2 iceRefractCoordR = texcoord.st - vec2(normalx * refractStrength, normaly * refractStrength) * (1.0f + abb * 2.0f);
		vec2 iceRefractCoordG = texcoord.st - vec2(normalx * refractStrength, normaly * refractStrength) * (1.0f + abb);
		vec2 iceRefractCoordB = texcoord.st - vec2(normalx * refractStrength, normaly * refractStrength) * (1.0f);
			
		vec3 refracted; 
		
					refracted.r	= texture2D(composite, iceRefractCoordR).r;
					refracted.g	= texture2D(composite, iceRefractCoordG).g;
					refracted.b	= texture2D(composite, iceRefractCoordB).b;
					
		color.rgb = refracted.rgb;
	
	}
*/	

#endif








#ifdef GLOSSY_REFLECTIONS

//color.rgb = texture2D(gnormal, texcoord.st).rgb;


const float glosslength = 0.020; //0.095
const float glosslength2 = 0.00575f;
const float g_distance = 50.0f;
const float fadefactor = 0.10f;
float gweight = 0.0f;
float gweight_add = 0.0f;
float gnormalpass = 0.0f;
float gdepthcheck = 0.0f;
const float gnoise = 0.0f;
float gnoise2 = 0.025f;
float gdepthdetect = 0.0f;
float gdepthr = texture2D(gdepth, texcoord.st).x;
float spreadx;
float spready;
float glossiness = 0.0f;


vec3 gloss = vec3(0.0f);
vec2 glosscoord = vec2(0.0f);
vec2 gcoord = texcoord.st;

float gfresnel = pow(distance(normalBiased, vec3(0.0, 0.0, 1.0))*1.5f, 1.0f);
	  
	  		float wetmask = clamp(texture2D(gaux2, texcoord.st).b * 4.0f - 2.5f, 0.0f, 1.0f) * wetx;
			float distmask = clamp(g_distance * fadefactor - linDepth * fadefactor, 0.0f, 1.0f);
			float g_spec = min(pow(specularity.r, 1.5f) * 1.4f + pow(specularity.g, 1.0f) * wetmask, 1.0f);
			float g_irr = texture2D(gaux3, texcoord.st).b * 0.95f;
			
			float totalspec = g_spec + g_irr;
			

		//for (int j = 0; j < 1; ++j) {
		
		//	glosscoord = vec2(0.0);
		//	spreadx = (0.5 - j*0.5) * 0.41;
		//	spreadx *= glossiness;
		//	spready = (0.5 - j*0.5) * 0.41;
		//	spready *= glossiness;
		
			for (int i = 0; i < 8; ++i) {
			
				if (linDepth > g_distance || totalspec < 0.001f) {
					break;
				}
				
				if (land == 0.0f) {
					break;
				}
						
				glosscoord += ((vec2(normalBiased.x + (noiseX2 * gnoise * normalBiased.x) * aspectRatio, normalBiased.y * 2.0f + (noiseX2 * gnoise * normalBiased.y * 2.0f)))*glosslength);
				//glosscoord += ((vec2(normalBiased.x + (noiseX2 * gnoise) * aspectRatio + spreadx, normalBiased.y * 2.0f + (noiseY2 * gnoise * 2.0f) + spready))*glosslength);
				
				//rotationMatrix = mat2( cos(gangle), sin(gangle),
									 // -sin(gangle), cos(gangle));
										
										
										
				//glosscoord *= rotationMatrix;
				
				gdepthdetect = texture2D(gdepth, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).x;

				gdepthcheck = (gdepthr - gdepthdetect);
				gdepthcheck = 1.0f - gdepthcheck;
				gdepthcheck = clamp(gdepthcheck * 480.0 - 479.0f, 0.0f, 1.0f);
				
				
				gnormalpass = clamp(distance(texture2D(gnormal, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).rgb, normal.rgb) * 2.5f - 0.5f, 0.0f, 2.0f);
				
				gweight = 1.0f * gnormalpass * gfresnel * gdepthcheck * (9 - i);
				
				gweight_add += gweight;


				
				gloss += max((texture2D(composite, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).rgb - 0.0f) * gweight, 0.0f);

			}		
		//}


			
			float finalweight =  (((gweight_add/105.0)*land) * (1.0f - iswater)) * 1.1;

			
			//color.rgb = color.rgb * (1.0f + wetmask * wetdark) - (wetdark * wetmask);
			
			gloss /= gweight_add * 1.1f;
			gloss = max(gloss, 0.0f);
			
			//reflect
			color.rgb = mix(color.rgb, gloss, clamp(finalweight * distmask * g_spec, 0.0f, 1.0f) * pow(dot(gloss, vec3(1.0f)), 0.2f));
			
			//metallic subtract
			//color.rgb *= (1.0f - g_irr * 0.75f);
			
			//metallic reflect
			color.rgb += gloss * albedo.rgb * finalweight * distmask * g_irr * pow(dot(gloss, vec3(1.0f)), 0.2f);
			//color.rgb = mix(color.rgb, gloss * color.rgb * 3.0f, finalweight * distmask * g_irr * pow(dot(gloss, vec3(1.0f)), 0.2f));
			//color.rgb = mix(color.rgb, gloss * albedo.rgb, clamp(finalweight * distmask * g_irr * albedo.rgb, 0.0f, 1.0f) * pow(dot(gloss, vec3(1.0f)), 0.2f));

			//radiosity
			//color.rgb *= 1.0f + (gloss * (1.0f - g_spec) * finalweight) * 0.5f;
			


			
#endif




#ifdef SCREEN_SPACE_RADIOSITY

//color.rgb = texture2D(gnormal, texcoord.st).rgb;

//vec3 normal = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;

const float radlength_r = 0.060;
const float radlength_r2 = 0.00575f;
float gweight_r = 0.0f;
float gnormalpass_r = 0.0f;
float gdepthcheck_r = 0.0f;
const float gnoise_r = 0.000f;
float gnoise_r2 = 0.125f;
float gdepthdetect_r = 0.0f;
float gweight_racc = 0.0f;


float sphere_2 = 0.0f;

vec3 rad = vec3(0.0f);
vec2 radcoord = vec2(0.0f);
vec2 gcoord = texcoord.st;

float gfresnel = pow(distance(normalBiased, vec3(0.0, 0.0, 1.0))*1.5f, 2.0f);
	  gfresnel = 4.0f;

			for (int i = 0; i < 5; ++i) {
			
				for (int k = 0; k < 6; ++k) {
				
					float rradius = pow(i, 1.5f);
			
					sphere_2 = (-1.0f + k + (mod(i, 2) - 1.0f))*rradius/2.0f;
						
					radcoord = ((vec2(normalBiased.x * rradius + normalBiased.y * (sphere_2) + (noiseX2 * gnoise_r) * aspectRatio, normalBiased.y * rradius * 2.0f + normalBiased.x * (sphere_2) * 2.0f + (noiseY2 * gnoise_r) * 2.0f))*radlength_r) / clamp(linDepth, 0.0f, 10.0f);
				
					gdepthcheck_r = distance(linDepth, getDepth(texcoord.st + radcoord)) / i;
					gdepthcheck_r = 1.0f - gdepthcheck_r;
					gdepthcheck_r = pow(clamp(gdepthcheck_r * 1.0 - 0.0f, 0.0f, 1.0f), 0.5f);
				
				
					gnormalpass_r = clamp(distance(texture2D(gnormal, gcoord.st + radcoord).rg, normal.rg) * 2.5f - 0.5f, 0.0f, 2.0f);
				
					gweight_r = 1.0f * gnormalpass_r * gfresnel * gdepthcheck_r  * pow((7 - i), 1.0f);
				
					gweight_racc += gweight_r;
					
					rad += max((texture2D(gcolor, gcoord.st + radcoord).rgb - 0.0f) * gweight_r, 0.0f);
				
				}
				
				if (land == 0.0f) {
					break;
				}
				

				


			}		


		
			
			float finalweight_r =  (((gweight_racc/25.0)*land) * (1.0f - iswater)) * 0.17 * RADIOSITY_AMOUNT;
			
			
			rad /= gweight_racc;
			rad = max(rad, 0.0f);
	
			//color.rgb += rad * finalweight_r * albedo * 0.15f;
			color.rgb *= 1.0 + ((rad - 0.000f) * finalweight_r * 1.25f);
			//color.rgb *= 1.0f - finalweight_r;
			
			//color.rgb *= max(1.0f - (finalweight_r * (1.0f - min(rad * 6.0f, 1.0f))), 0.0f);

#endif





#ifdef GODRAYS

	float GR = addGodRays(0.0f, Texcoord2, noiseX3, noiseX4, noiseY4, noiseX2, noiseY2, noiseX5, noiseY5, noiseX6, noiseY6)/2.0;

	float GRr = 1.0 - texture2D(gaux2, texcoord.st).r;
	
	//GR = mix(GR, 0.0f, rainx);
	
	/*
	float GRs  = 1.0 - texture2D(gaux1, vec2(0.55, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.55, 0.45)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.45)).g;

		  GRs /= 3.0;
	*/
	
	vec3 sunrise_sun;
	 sunrise_sun.r = 1.0 * TimeSunrise;
	 sunrise_sun.g = 0.629 * TimeSunrise;
	 sunrise_sun.b = 0.416 * TimeSunrise;
	
	vec3 noon_sun;
	 noon_sun.r = 1.0 * TimeNoon;
	 noon_sun.g = 1.0 * TimeNoon;
	 noon_sun.b = 0.98 * TimeNoon;
	
	vec3 sunset_sun;
	 sunset_sun.r = 0.99 * TimeSunset;
	 sunset_sun.g = 0.839 * TimeSunset;
	 sunset_sun.b = 0.666 * TimeSunset;
	
	vec3 midnight_sun;
	 midnight_sun.r = 0.45 * TimeMidnight * 0.20f;
	 midnight_sun.g = 0.70 * TimeMidnight * 0.20f;
	 midnight_sun.b = 1.00 * TimeMidnight * 0.20f;
	 
	vec3 rain_sun_day;
	 rain_sun_day.r = 1.0f * (1.0f - TimeMidnight) * 0.1f; 
	 rain_sun_day.g = 1.0f * (1.0f - TimeMidnight) * 0.1f;
	 rain_sun_day.b = 1.0f * (1.0f - TimeMidnight) * 0.1f;	
	 
	vec3 rain_sun_night;
	 rain_sun_night.r = 1.0f * (TimeMidnight) * 0.0f;
	 rain_sun_night.g = 1.0f * (TimeMidnight) * 0.0f;
	 rain_sun_night.b = 1.0f * (TimeMidnight) * 0.0f;
	
	vec3 sunlight;
	 sunlight = mix(sunrise_sun + noon_sun + sunset_sun + midnight_sun, rain_sun_day + rain_sun_night, rainx);

	
	

	
	
	GR = pow(GR, 1.0f)*2.5f;
	
	color.r += pow(GR*sunlight.r, 1.0f);
	color.g += pow(GR*sunlight.g, 1.0f);
	color.b += pow(GR*sunlight.b, 1.0f);
	
	
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
	
	float radius = 0.002f*GLARE_RANGE;
	const float radiusv = 0.002f;
	const float bloomintensity = 0.1f*GLARE_AMOUNT;
	
	const float glarex = 2.0f;
	const float glaresub = 0.0f;
	
	float bloomnoise = noiseX2*0.0f;
	float bloomnoisey = noiseY2*0.0f;
	

	vec4 clr = vec4(0.0f);
	
	//clr += texture2D(composite, texcoord.st);
	
	for (int i = 0; i < 1; ++i) {
	//horizontal (70 taps)

	clr +=  max(texture2D(composite, texcoord.st + (vec2(5.0f+bloomnoise,5.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(4.0f+bloomnoise,4.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(3.0f+bloomnoise,3.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(2.0f+bloomnoise,2.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(1.0f+bloomnoise,1.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
		//clr += texture2D(composite, texcoord.st + (vec2(0.0f,0.0f))*radius)*6.0f;
		
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-1.0f+bloomnoise,1.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-2.0f+bloomnoise,2.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-3.0f+bloomnoise,3.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-4.0f+bloomnoise,4.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-5.0f+bloomnoise,5.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;

	//vertical

	clr +=  max(texture2D(composite, texcoord.st + (vec2(5.0+bloomnoise,-5.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(4.0+bloomnoise,-4.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(3.0+bloomnoise,-3.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(2.0+bloomnoise,-2.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(1.0+bloomnoise,-1.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-5.0+bloomnoise,-5.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-4.0+bloomnoise,-4.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-3.0+bloomnoise,-3.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-2.0+bloomnoise,-2.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(composite, texcoord.st + (vec2(-1.0+bloomnoise,-1.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
	radius *= 3.0;
	
	clr = (clr/10.0f)/2.0f;
	}
	
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
	
	float radius2 = 0.006f*GLARE_RANGE;
	const float radius2v = 0.002f;
	const float bloomintensity2 = 0.08f*GLARE_AMOUNT;
	
	const float glarex2 = 2.0f;
	const float glaresub2 = 0.0f;
	
	float bloomnoise2 = noiseX4*0.0f;	
	float bloomnoise2y = noiseY4*0.0f;	

	vec4 clr2 = vec4(0.0f);
	
	//clr2 += texture2D(composite, texcoord.st);
	
	//horizontal (70 taps)
	
	for (int i = 0; i < 1; ++i) {

	clr2 += max(texture2D(composite, texcoord.st + (vec2(5.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(4.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(3.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(2.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(1.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
		//clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0f,0.0f))*radius2)*6.0f;
		
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-1.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-2.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-3.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-4.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(-5.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;

	//vertical

	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-5.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-4.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-3.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-2.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,-1.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,5.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,4.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,3.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,2.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(composite, texcoord.st + (vec2(0.0+bloomnoise2,1.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
	radius2 /= 3.0;
	
	clr2 = (clr2/10.0f)/2.0f;
	}
	
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
			//lPos = clamp(lPos, vec2(0.001f), vec2(0.999f));
			

			
			float sunmask = 0.0f;
			float sunstep = -4.5f;
			float masksize = 0.004f;
					

					sunmask += 1.0f - texture2D(gaux1, lPos).b;
					
					if (lPos.x > 1.0f || lPos.x < 0.0f || lPos.y > 1.0f || lPos.y < 0.0f) {
							sunmask = 0.0f;
					}

					sunmask *= LENS_POWER * (1.0f - TimeMidnight);
					sunmask *= 1.0 - rainx;
			
			//Detect if sun is on edge of screen
				float edgemaskx = clamp(distance(lPos.x, 0.5f)*8.0f - 3.0f, 0.0f, 1.0f);
				float edgemasky = clamp(distance(lPos.y, 0.5f)*8.0f - 3.0f, 0.0f, 1.0f);
			
						
						
			////Darken colors if the sun is visible
				float centermask = 1.0 - clamp(distance(lPos.xy, vec2(0.5f, 0.5f))*2.0, 0.0, 1.0);
						centermask = pow(centermask, 1.0f);
						centermask *= sunmask;
			
				color.r *= (1.0 - centermask * (1.0f - TimeMidnight));
				color.g *= (1.0 - centermask * (1.0f - TimeMidnight));
				color.b *= (1.0 - centermask * (1.0f - TimeMidnight));
			
			
			//Adjust global flare settings
				const float flaremultR = 0.8f;
				const float flaremultG = 1.0f;
				const float flaremultB = 1.5f;
			
				float flarescale = 1.0f;
				const float flarescaleconst = 1.0f;
			
			
			//Flare gets bigger at center of screen
			
				flarescale *= (1.0 - centermask);
			

			//Center white flare
			vec2 flare1scale = vec2(1.7f*flarescale, 1.7f*flarescale);
			float flare1pow = 12.0f;
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
			  vec2 flare1Bscale = vec2(0.5f*flarescale, 0.5f*flarescale);
			  float flare1Bpow = 6.0f;
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
			  vec2 flare3scale = vec2(2.0f*flarescale, 2.0f*flarescale);
			  float flare3pow = 0.7f;
			  float flare3fill = 10.0f;
			  float flare3offset = -0.5f;
			vec2 flare3pos = vec2(  ((1.0 - lPos.x)*(flare3offset + 1.0) - (flare3offset*0.5))  *aspectRatio*flare3scale.x,  ((1.0 - lPos.y)*(flare3offset + 1.0) - (flare3offset*0.5))  *flare3scale.y);
			
			
			float flare3 = distance(flare3pos, vec2(texcoord.s*aspectRatio*flare3scale.x, texcoord.t*flare3scale.y));
				  flare3 = 0.5 - flare3;
				  flare3 = clamp(flare3*flare3fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare3 = sin(flare3*1.57075);
				  flare3 *= sunmask;
				  flare3 = pow(flare3, 1.1f);
				  
				  flare3 *= flare3pow;			
				  
				  
				  //subtract from blue flare
				  vec2 flare3Bscale = vec2(1.4f*flarescale, 1.4f*flarescale);
				  float flare3Bpow = 1.0f;
				  float flare3Bfill = 2.0f;
				  float flare3Boffset = -0.65f;
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
			  vec2 flare3Cscale = vec2(3.2f*flarescale, 3.2f*flarescale);
			  float flare3Cpow = 1.4f;
			  float flare3Cfill = 10.0f;
			  float flare3Coffset = -0.0f;
			vec2 flare3Cpos = vec2(  ((1.0 - lPos.x)*(flare3Coffset + 1.0) - (flare3Coffset*0.5))  *aspectRatio*flare3Cscale.x,  ((1.0 - lPos.y)*(flare3Coffset + 1.0) - (flare3Coffset*0.5))  *flare3Cscale.y);
			
			
			float flare3C = distance(flare3Cpos, vec2(texcoord.s*aspectRatio*flare3Cscale.x, texcoord.t*flare3Cscale.y));
				  flare3C = 0.5 - flare3C;
				  flare3C = clamp(flare3C*flare3Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare3C = sin(flare3C*1.57075);
				  
				  flare3C = pow(flare3C, 1.1f);
				  
				  flare3C *= flare3Cpow;			
				  
				  
				  //subtract from blue flare
				  vec2 flare3Dscale = vec2(2.1f*flarescale, 2.1f*flarescale);
				  float flare3Dpow = 2.7f;
				  float flare3Dfill = 1.4f;
				  float flare3Doffset = -0.05f;
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
			  vec2 flare4scale = vec2(4.5f*flarescale, 4.5f*flarescale);
			  float flare4pow = 0.3f;
			  float flare4fill = 3.0f;
			  float flare4offset = -0.1f;
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
			  vec2 flare4Bscale = vec2(7.5f*flarescale, 7.5f*flarescale);
			  float flare4Bpow = 0.4f;
			  float flare4Bfill = 2.0f;
			  float flare4Boffset = 0.0f;
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
			  vec2 flare4Cscale = vec2(37.5f*flarescale, 37.5f*flarescale);
			  float flare4Cpow = 2.0f;
			  float flare4Cfill = 2.0f;
			  float flare4Coffset = -0.3f;
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
			  vec2 flare4Dscale = vec2(67.5f*flarescale, 67.5f*flarescale);
			  float flare4Dpow = 1.0f;
			  float flare4Dfill = 2.0f;
			  float flare4Doffset = -0.35f;
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
			  vec2 flare4Escale = vec2(60.5f*flarescale, 60.5f*flarescale);
			  float flare4Epow = 1.0f;
			  float flare4Efill = 3.0f;
			  float flare4Eoffset = -0.3393f;
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
			  vec2 flare4Fscale = vec2(20.5f*flarescale, 20.5f*flarescale);
			  float flare4Fpow = 3.0f;
			  float flare4Ffill = 3.0f;
			  float flare4Foffset = -0.4713f;
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
			  vec2 flare5scale = vec2(3.2f*flarescale , 3.2f*flarescale );
			  float flare5pow = 13.4f;
			  float flare5fill = 1.0f;
			  float flare5offset = -2.0f;
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
			  vec2 flare5Bscale = vec2(0.5f*flarescale , 3.5f*flarescale );
			  float flare5Bpow = 1.4f;
			  float flare5Bfill = 2.0f;
			  float flare5Boffset = -1.9f;
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
			  vec2 flare6scale = vec2(1.2f*flarescale, 1.2f*flarescale);
			  float flare6pow = 0.2f;
			  float flare6fill = 5.0f;
			  float flare6offset = -1.9f;
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
			  vec2 flare6Bscale = vec2(1.1f*flarescale, 1.1f*flarescale);
			  float flare6Bpow = 0.2f;
			  float flare6Bfill = 5.0f;
			  float flare6Boffset = -1.9f;
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
			  vec2 flare6Cscale = vec2(0.9f*flarescale, 0.9f*flarescale);
			  float flare6Cpow = 0.2f;
			  float flare6Cfill = 5.0f;
			  float flare6Coffset = -1.9f;
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

			  vec2 flare7scale = vec2(0.4f*flarescale, 0.4f*flarescale);
			  float flare7pow = 0.2f;
			  float flare7fill = 10.0f;
			  float flare7offset = 2.6f;
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

			  vec2 flare7Dscale = vec2(0.39f*flarescale, 0.39f*flarescale);
			  float flare7Dpow = 0.1f;
			  float flare7Dfill = 10.0f;
			  float flare7Doffset = 2.6f;
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

			  vec2 flare7Bscale = vec2(0.2f*flarescale, 0.2f*flarescale);
			  float flare7Bpow = 0.1f;
			  float flare7Bfill = 2.0f;
			  float flare7Boffset = 2.9f;
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
			  vec2 flare8scale = vec2(0.3f*flarescale, 40.5f*flarescale);
			  float flare8pow = 0.5f;
			  float flare8fill = 12.0f;
			  float flare8offset = 1.0f;
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
			  vec2 flare9scale = vec2(0.2f*flarescale, 5.5f*flarescale);
			  float flare9pow = 1.9f;
			  float flare9fill = 2.0f;
			  vec2 flare9offset = vec2(1.0f, 0.0f);
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
			  vec2 flare10scale = vec2(6.0f*flarescale, 6.0f*flarescale);
			  float flare10pow = 1.9f;
			  float flare10fill = 1.1f;
			  float flare10offset = -0.7f;
			vec2 flare10pos = vec2(  ((1.0 - lPos.x)*(flare10offset + 1.0) - (flare10offset*0.5))  *aspectRatio*flare10scale.x,  ((1.0 - lPos.y)*(flare10offset + 1.0) - (flare10offset*0.5))  *flare10scale.y);
			
			
			float flare10 = distance(flare10pos, vec2(texcoord.s*aspectRatio*flare10scale.x, texcoord.t*flare10scale.y));
				  flare10 = 0.5 - flare10;
				  flare10 = clamp(flare10*flare10fill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare10 = sin(flare10*1.57075);
				  flare10 *= sunmask;
				  flare10 = pow(flare10, 1.1f);
				  
				  flare10 *= flare10pow;			
				  
				  
				  //subtract
				  vec2 flare10Bscale = vec2(5.1f*flarescale, 5.1f*flarescale);
				  float flare10Bpow = 1.5f;
				  float flare10Bfill = 1.0f;
				  float flare10Boffset = -0.77f;
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
			  vec2 flare10Cscale = vec2(6.0f*flarescale, 6.0f*flarescale);
			  float flare10Cpow = 1.9f;
			  float flare10Cfill = 1.1f;
			  float flare10Coffset = -0.6f;
			vec2 flare10Cpos = vec2(  ((1.0 - lPos.x)*(flare10Coffset + 1.0) - (flare10Coffset*0.5))  *aspectRatio*flare10Cscale.x,  ((1.0 - lPos.y)*(flare10Coffset + 1.0) - (flare10Coffset*0.5))  *flare10Cscale.y);
			
			
			float flare10C = distance(flare10Cpos, vec2(texcoord.s*aspectRatio*flare10Cscale.x, texcoord.t*flare10Cscale.y));
				  flare10C = 0.5 - flare10C;
				  flare10C = clamp(flare10C*flare10Cfill, 0.0, 1.0) * clamp(-sP.z, 0.0, 1.0);
				  flare10C = sin(flare10C*1.57075);
				  flare10C *= sunmask;
				  flare10C = pow(flare10C, 1.1f);
				  
				  flare10C *= flare10Cpow;			
				  
				  
				  //subtract
				  vec2 flare10Dscale = vec2(5.1f*flarescale, 5.1f*flarescale);
				  float flare10Dpow = 1.5f;
				  float flare10Dfill = 1.0f;
				  float flare10Doffset = -0.67f;
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

			  vec2 flare11scale = vec2(1.5f*flarescale, 1.5f*flarescale);
			  float flare11pow = 1.1f;
			  float flare11fill = 2.0f;
			  float flare11offset = -0.523f;
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

			  vec2 flare12scale = vec2(2.5f*flarescale, 2.5f*flarescale);
			  float flare12pow = 0.5f;
			  float flare12fill = 2.0f;
			  float flare12offset = -0.323f;
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

			  vec2 flare13scale = vec2(1.0f*flarescale, 1.0f*flarescale);
			  float flare13pow = 0.5f;
			  float flare13fill = 2.0f;
			  float flare13offset = +0.138f;
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


#endif

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
	//color = max(((color*1.10f) - 0.06f), 0.0f);

	//color.r = pow(color.r, GAMMA);
	//color.g = pow(color.g, GAMMA);
	//color.b = pow(color.b, GAMMA);

	

	
//color *= 1.1f;

#ifdef VINTAGE

	color.r = clamp(color.r, 0.04, 1.0);

	color.b = clamp(color.b, 0.06, 0.89);
	

#endif

#ifdef LOWLIGHT_EYE
	
	vec3 rodcolor = mix(vec3(0.1f, 0.25f, 1.0f), vec3(1.0f), 0.5f);

	color.rgb = mix(color.rgb, vec3(dot(color.rgb, vec3(1.0f))) * rodcolor, clamp(1.0f - dot(color.rgb, vec3(1.0f)) * 8.0f, 0.0f, 1.0f) * 0.5f);

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



	color.rgb *= 0.8f;
	 float TonemapOversat = 4.4f;
	 float TonemapCurve   = 0.17f;
	 
	#ifdef TONEMAP_FILMIC

		TonemapOversat = 5.0f;
		TonemapCurve   = 0.199f;
		//#define COLOR_BOOST	0.09
		
		#ifdef TONEMAP_COLOR_FILTER
			color.r -= 0.01f;
			color.g -= 0.005f;
			color.b = color.b * 0.9f + 0.00f;
			
			
			
		#endif

	
	#endif
	
	//Tonemap
	color.rgb = (color.rgb * (1.0 + color.rgb/TonemapOversat))/(color.rgb + TonemapCurve);
	
	color = color*(1.0f + DARKMULT) - DARKMULT;
	
	color = clamp(color, 0.0f, 1.0f);

#endif
/*
float exposureb = 3.0f;

color.rgb = (1.0f - exp( -color.rgb * exposureb ));
*/

color.rgb *= 1.25f;

//color.rgb = pow(color.rgb, vec3(1/1.5));


	//Color boosting
	color.r = (color.r)*(COLOR_BOOST + 1.0f) + (color.g + color.b)*(-COLOR_BOOST);
	color.g = (color.g)*(COLOR_BOOST + 1.0f) + (color.r + color.b)*(-COLOR_BOOST);
	color.b = (color.b)*(COLOR_BOOST + 1.0f) + (color.r + color.g)*(-COLOR_BOOST);
	
	color.r = pow(color.r, GAMMA);
	color.g = pow(color.g, GAMMA);
	color.b = pow(color.b, GAMMA);
	
	gl_FragColor = color;
	
// End of Main. -----------------
}
