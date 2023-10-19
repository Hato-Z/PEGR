#version 120





//to increase shadow draw distance, edit SHADOWDISTANCE and SHADOWHPL below. Both should be equal. Needs decimal point.
//disabling is done by adding "//" to the beginning of a line.





//ADJUSTABLE VARIABLES

#define BLURFACTOR 3.5
#define SHADOW_DARKNESS 1.650   // 1.0 Is defualt darkness. 2.0 is black shadows. 0.0 is no shadows.
#define SHADOWDISTANCE 80.0 
#define SHADOW_CLAMP 0.4
#define SHADOW_RES 1024

/* SHADOWRES:1024 */
/* SHADOWHPL:80.0 */

  #define SSAO
  #define SSAO_LUMINANCE 0.0				// At what luminance will SSAO's shadows become highlights.
  #define SSAO_STRENGTH 1.75               // Too much strength causes white highlights on extruding edges and behind objects
  #define SSAO_LOOP 1						// Integer affecting samples that are taken to calculate SSAO. Higher values mean more accurate shadowing but bigger performance impact
  #define SSAO_NOISE true					// Randomize SSAO sample gathering. With noise enabled and SSAO_LOOP set to 1, you will see higher performance at the cost of fuzzy dots in shaded areas.
  #define SSAO_NOISE_AMP 0.0					// Multiplier of noise. Higher values mean SSAO takes random samples from a larger radius. Big performance hit at higher values.
  #define SSAO_MAX_DEPTH 0.9				// View distance of SSAO
  #define SSAO_SAMPLE_DELTA 0.4			// Radius of SSAO shadows. Higher values cause more performance hit.
  #define CORRECTSHADOWCOLORS				// Colors sunlight and ambient light correctly according to real-life. 
  #define SHADOWOFFSET 0.0				// Shadow offset multiplier. Values that are too low will cause artefacts.
  //#define FXAA							// FXAA shader. Broken, but you can give it a try if you want.
  #define GODRAYS
  #define GODRAYS_EXPOSURE 0.10
  #define GODRAYS_SAMPLES 6
  #define GODRAYS_DECAY 0.95
  #define GODRAYS_DENSITY 0.65

  #define SKY_LIGHTING
  #define SKY_LIGHTING_SPREAD 2.63f
  #define SKY_LIGHTING_MIN_DARKNESS 0.12f
  
  #define SKY_DESATURATION 0.1
  
  //#define CLAY_RENDER

//#define PRESERVE_COLOR_RANGE


//END OF ADJUSTABLE VARIABLES






uniform sampler2D gcolor;
uniform sampler2D gdepth;
uniform sampler2D gnormal;
uniform sampler2D composite;
uniform sampler2D shadow;
uniform sampler2D gaux1;
uniform sampler2D gaux2;
uniform sampler2D gaux3;

varying vec4 texcoord;
varying vec4 lmcoord;
varying vec3 lightVector;

uniform int worldTime;

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferModelViewInverse;
uniform mat4 shadowProjection;
uniform mat4 shadowModelView;
uniform vec3 sunPosition;

//attribute vec4 mc_Entity;

uniform float near;
uniform float far;
uniform float viewWidth;
uniform float viewHeight;
uniform float rainStrength;
uniform float wetness;
uniform float aspectRatio;



// Standard depth function.
float getDepth(vec2 coord) {
    return 2.0f * near * far / (far + near - (2.0f * texture2D(gdepth, coord).x - 1.0f) * (far - near));
}

//Auxilliary variables
float	land 			 = texture2D(gaux1, texcoord.st).b;
//float	noblur 			 = texture2D(gaux1, texcoord.st).r;
vec3	sunPos			 = sunPosition;
vec2 	Texcoord2		 = texcoord.st;
float 	iswater			 = texture2D(gaux1, texcoord.st).g;
vec3 	normal			 = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
float 	translucent		 = texture2D(gaux1, texcoord.st).r;
float   texshading		 = 1.0f;
float   isice 			 = texture2D(gaux3, texcoord.st).b;
float	specularityDry	 = texture2D(gaux3, texcoord.st).r;
float 	specularityWet	 = texture2D(gaux3, texcoord.st).g;

//Crossfading conditionals

float rainx = clamp(rainStrength, 0.0f, 1.0f);
float wetx  = clamp(wetness, 0.0f, 1.0f);
float landx = land;

//Lightmaps

float sky_lightmap = texture2D(gaux2, texcoord.st).b;
float torch_lightmap = texture2D(gaux2, texcoord.st).r;
float lightning_lightmap = texture2D(gaux2, texcoord.st).g;


//Calculate Time of Day

	float timefract = worldTime;
	float timePow = 2.0f;

	float TimeSunrise  = ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0) + (1.0 - (clamp(timefract, 0.0, 6000.0)/6000.0));
		  
	float TimeNoon     = ((clamp(timefract, 0.0, 6000.0)) / 6000.0) - ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0);
	  
	float TimeSunset   = ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0) - ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0);
		  
	float TimeMidnight = ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0) - ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0);




#ifdef SSAO

// Alternate projected depth (used by SSAO, probably AA too)
float getProDepth(vec2 coord) {
	float depth = texture2D(gdepth, coord).x;
	return ( 2.0f * near ) / ( far + near - depth * ( far - near ) );
}

float znear = near; //Z-near
float zfar = far; //Z-far

float diffarea = 0.6f; //self-shadowing reduction
float gdisplace = 0.30f; //gauss bell center

//bool noise = SSAO_NOISE; //use noise instead of pattern for sample dithering?
bool onlyAO = false; //use only ambient occlusion pass?

vec2 texCoord = texcoord.st;


vec2 rand(vec2 coord) { //generating noise/pattern texture for dithering
  const float width = 1.0f;
  const float height = 1.0f;
  float noiseX = ((fract(1.0f-coord.s*(width/2.0f))*0.25f)+(fract(coord.t*(height/2.0f))*0.75f))*2.0f-1.0f;
  float noiseY = ((fract(1.0f-coord.s*(width/2.0f))*0.75f)+(fract(coord.t*(height/2.0f))*0.25f))*2.0f-1.0f;

  //generate SSAO noise
  noiseX = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f))) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  noiseY = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f)*2.0f)) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  
  return vec2(noiseX,noiseY)*0.002f*SSAO_NOISE_AMP;
}


float compareDepths(in float depth1, in float depth2, in int zfar) {  
  float garea = 8.5f; //gauss bell width    
  float diff = (depth1 - depth2) * 100.0f; //depth difference (0-100)
  //reduce left bell width to avoid self-shadowing 
  
  if (diff < gdisplace) {
    garea = diffarea;
  } else {
    zfar = 1;
  }


  float gauss = pow(2.7182f,-2.0f*(diff-gdisplace)*(diff-gdisplace)/(garea*garea));
  return gauss;
} 

float calAO(float depth, float dw, float dh) {  
  float temp = 0.0f;
  float temp2 = 0.0f;
  dw *= 2.0f;
  dh *= 2.0f;
  float coordw = texCoord.x + dw/(depth*0.2f + 0.1f);
  float coordh = texCoord.y + dh/(depth*0.2f + 0.1f);
  float coordw2 = texCoord.x - dw/(depth*0.2f + 0.1f);
  float coordh2 = texCoord.y - dh/(depth*0.2f + 0.1f);

  if (coordw  < 1.0f && coordw  > 0.0f && coordh < 1.0f && coordh  > 0.0f){
    vec2 coord = vec2(coordw , coordh);
    vec2 coord2 = vec2(coordw2, coordh2);
    int zfar = 0;
    temp = compareDepths(depth, getProDepth(coord),zfar);

    //DEPTH EXTRAPOLATION:
    //if (zfar > 0){
    //  temp2 = compareDepths(getProDepth(coord2),depth,zfar);
    //  temp += (1.0f-temp)*temp2; 
    //}
  }

  return temp;  
}  



float getSSAOFactor() {

  float incx = 1.0f / viewWidth * SSAO_SAMPLE_DELTA;
  float incy = 1.0f / viewHeight * SSAO_SAMPLE_DELTA;
  
  
	vec2 noise1 = rand(texCoord)*20.0f; 
	
	/*
	vec2 noise2 = rand(texCoord + vec2(incx, incy)*10); 
	vec2 noise3 = rand(texCoord + vec2(incx, -incy)*10); 
	vec2 noise4 = rand(texCoord + vec2(-incx, incy)*10); 
	vec2 noise5 = rand(texCoord + vec2(-incx, -incy)*10); 
	*/
	
	
	float depth = getProDepth(texCoord);
  if (depth > SSAO_MAX_DEPTH) {
    return 1.0f;
  }
  float cdepth = texture2D(gdepth,texCoord).g;
	
	float ao;
	float s;
	

  float pw = incx;
  float ph = incy;
  float aoMult = SSAO_STRENGTH;
  int aaLoop = SSAO_LOOP;
  float aaDiff = (1.0f + 2.0f / 1.0f); // 1.0 is samples

    float npw  = (pw + 0.05f * noise1.x) / cdepth;
    float nph  = (ph + 0.05f * noise1.y) / cdepth;
	

	float npw2  = (pw*2.0f + 0.05f * noise1.x) / cdepth;
    float nph2  = (ph*2.0f + 0.05f * noise1.y) / cdepth;
	
	float npw3  = (pw*3.0f + 0.05f * noise1.x) / cdepth;
    float nph3  = (ph*3.0f + 0.05f * noise1.y) / cdepth;
	
	float npw4  = (pw*4.0f + 0.05f * noise1.x) / cdepth;
    float nph4  = (ph*4.0f + 0.05f * noise1.y) / cdepth;

    ao += calAO(depth, npw, nph) * aoMult;
    ao += calAO(depth, npw, -nph) * aoMult;
    ao += calAO(depth, -npw, nph) * aoMult;
    ao += calAO(depth, -npw, -nph) * aoMult;
	
	ao += calAO(depth, npw2, nph2) * aoMult/1.5f;
    ao += calAO(depth, npw2, -nph2) * aoMult/1.5f;
    ao += calAO(depth, -npw2, nph2) * aoMult/1.5f;
    ao += calAO(depth, -npw2, -nph2) * aoMult/1.5f;
	
	ao += calAO(depth, npw3, nph3) * aoMult/2.0f;
    ao += calAO(depth, npw3, -nph3) * aoMult/2.0f;
    ao += calAO(depth, -npw3, nph3) * aoMult/2.0f;
    ao += calAO(depth, -npw3, -nph3) * aoMult/2.0f;
	
	
	/*
	ao += calAO(depth, npw4, nph4) * aoMult/2.5f;
    ao += calAO(depth, npw4, -nph4) * aoMult/2.5f;
    ao += calAO(depth, -npw4, nph4) * aoMult/2.5f;
    ao += calAO(depth, -npw4, -nph4) * aoMult/2.5f;
	
	
	 ao += calAO(depth, 2.0*npw2, 2.0*nph2) * aoMult/2.0;
    ao += calAO(depth, 2.0*npw2, -2.0*nph2) * aoMult/2.0;
    ao += calAO(depth, -2.0*npw2, 2.0*nph2) * aoMult/2.0;
    ao += calAO(depth, -2.0*npw2, -2.0*nph2) * aoMult/2.0;
	
	 ao += calAO(depth, 3.0*npw3, 3.0*nph3) * aoMult/3.0;
    ao += calAO(depth, 3.0*npw3, -3.0*nph3) * aoMult/3.0;
    ao += calAO(depth, -3.0*npw3, 3.0*nph3) * aoMult/3.0;
    ao += calAO(depth, -3.0*npw3, -3.0*nph3) * aoMult/3.0;
	
	 ao += calAO(depth, 4.0*npw4, 4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, 4.0*npw4, -4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, -4.0*npw4, 4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, -4.0*npw4, -4.0*nph4) * aoMult/4.0;
	*/
	
	ao /= 16.0f;
	ao = 1.0f-ao;	
  ao = clamp(ao, 0.0f, 0.5f) * 2.0f;
	
  return ao;
}

#endif


#ifdef GODRAYS



	float addGodRays(in float nc, in vec2 tx, in float noise, in float noise2, in float noise3, in float noise4, in float noise5) {
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
			
				tx -= delta;
				float sample = 0.0f;

					sample = 1.0f - texture2D(gaux1, tx + vec2(noise*delta.x, noise*delta.y)).b;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise2*delta.x, noise2*delta.y)).b;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise3*delta.x, noise3*delta.y)).b;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise4*delta.x, noise4*delta.y)).b;
					sample += 1.0f - texture2D(gaux1, tx + vec2(noise5*delta.x, noise5*delta.y)).b;
				sample *= decay;

					colorGD += sample;
					decay *= GODRAYS_DECAY;
			}
			
			float bubble = distance(vec2(delta.x*aspectRatio, delta.y), vec2(0.0f, 0.0f))*4.0f;
				  bubble = clamp(bubble, 0.0f, 1.0f);
				  bubble = 1.0f - bubble;
				  
			return (nc + GODRAYS_EXPOSURE * (colorGD*bubble))*GDTimeMult;
        
	}
#endif 









///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void main() {

land 			 = texture2D(gaux1, texcoord.st).b;
landx			 = land;
//noblur 		 = texture2D(gaux1, texcoord.st).r;
iswater			 = texture2D(gaux1, texcoord.st).g;
normal         	 = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
translucent		 = texture2D(gaux1, texcoord.st).r;
isice 			 = texture2D(gaux3, texcoord.st).b;
specularityDry   = texture2D(gaux3, texcoord.st).r;
specularityWet   = texture2D(gaux3, texcoord.st).g;

//Lightmaps

sky_lightmap = texture2D(gaux2, texcoord.st).b;
torch_lightmap = texture2D(gaux2, texcoord.st).r;
lightning_lightmap = texture2D(gaux2, texcoord.st).g;

//Curve times
//Curve times
//Curve times
		  TimeSunrise  = pow(TimeSunrise, timePow);
		  TimeNoon     = pow(TimeNoon, 1.0f/timePow);
		  TimeSunset   = pow(TimeSunset, timePow);
		  TimeMidnight = pow(TimeMidnight, 1.0f/timePow);

float noiseamp = 0.3f;
					
						float width2 = 1.0f;
						float height2 = 1.0f;
						float noiseX2 = ((fract(1.0f-Texcoord2.s*(width2/2.0f))*0.25f)+(fract(Texcoord2.t*(height2/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY2 = ((fract(1.0f-Texcoord2.s*(width2/2.0f))*0.75f)+(fract(Texcoord2.t*(height2/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX2 = clamp(fract(sin(dot(Texcoord2 ,vec2(12.9898f,78.233f))) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY2 = clamp(fract(sin(dot(Texcoord2 ,vec2(12.9898f,78.233f)*2.0f)) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX2 *= (0.0005f*noiseamp);
						noiseY2 *= (0.0005f*noiseamp);
						
						float width3 = 2.0f;
						float height3 = 2.0f;
						float noiseX3 = ((fract(1.0f-Texcoord2.s*(width3/2.0f))*0.25f)+(fract(Texcoord2.t*(height3/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY3 = ((fract(1.0f-Texcoord2.s*(width3/2.0f))*0.75f)+(fract(Texcoord2.t*(height3/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX3 = clamp(fract(sin(dot(Texcoord2 ,vec2(18.9898f,28.633f))) * 4378.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY3 = clamp(fract(sin(dot(Texcoord2 ,vec2(11.9898f,59.233f)*2.0f)) * 3758.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX3 *= (0.0005f*noiseamp);
						noiseY3 *= (0.0005f*noiseamp);
						
						float width4 = 3.0f;
						float height4 = 3.0f;
						float noiseX4 = ((fract(1.0f-Texcoord2.s*(width4/2.0f))*0.25f)+(fract(Texcoord2.t*(height4/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY4 = ((fract(1.0f-Texcoord2.s*(width4/2.0f))*0.75f)+(fract(Texcoord2.t*(height4/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX4 = clamp(fract(sin(dot(Texcoord2 ,vec2(16.9898f,38.633f))) * 41178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY4 = clamp(fract(sin(dot(Texcoord2 ,vec2(21.9898f,66.233f)*2.0f)) * 9758.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX4 *= (0.0005f*noiseamp);
						noiseY4 *= (0.0005f*noiseamp);
						
						float width5 = 4.0f;
						float height5 = 4.0f;
						float noiseX5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.25f)+(fract(Texcoord2.t*(height5/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.75f)+(fract(Texcoord2.t*(height5/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX5 = clamp(fract(sin(dot(Texcoord2 ,vec2(11.9898f,68.633f))) * 21178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY5 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898f,71.233f)*2.0f)) * 6958.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX5 *= (0.0005f*noiseamp);
						noiseY5 *= (0.0005f*noiseamp);
						

//



	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * texture2D(gdepth, texcoord.st).y - 1.0f, 1.0f);
	fragposition /= fragposition.w;
	
	#ifdef SHADOWDISTANCE
	float drawdistance = SHADOWDISTANCE;
	float drawdistancesquared = pow(drawdistance, 2.0f);
	#endif
	
	float distance = sqrt(fragposition.x * fragposition.x + fragposition.y * fragposition.y + fragposition.z * fragposition.z);

	float shading = 1.0f;
	float shadingsharp = 1.0f;
	float shadingao = 1.0f;
	
	
	vec4 worldposition = vec4(0.0);
	vec4 worldpositionraw = vec4(0.0);
			
	worldposition = gbufferModelViewInverse * fragposition;	
	
	float xzDistanceSquared = worldposition.x * worldposition.x + worldposition.z * worldposition.z;
	float yDistanceSquared  = worldposition.y * worldposition.y;
	
	worldpositionraw = worldposition;
	
			worldposition = shadowModelView * worldposition;
			float comparedepth = -worldposition.z;
			worldposition = shadowProjection * worldposition;
			worldposition /= worldposition.w;
			
			worldposition.st = worldposition.st * 0.5f + 0.5f;
			
			
			
			
			////////////////////////////////////WAVES////////////////////////////
			////////////////////////////////////WAVES////////////////////////////
			////////////////////////////////////WAVES////////////////////////////
float wsize = 9.0f*3.0;
float wspeed = 0.3f;

float rs0 = abs(sin((worldTime*wspeed/5.0) + (worldposition.s*wsize) * 20.0 + (worldposition.z*4.0))+0.2);
float rs1 = abs(sin((worldTime*wspeed/7.0) + (worldposition.t*wsize) * 27.0) + 0.5);
float rs2 = abs(sin((worldTime*wspeed/2.0) + (worldposition.t*wsize) * 60.0 - sin(worldposition.s*wsize) * 13.0)+0.4);
float rs3 = abs(sin((worldTime*wspeed/1.0) - (worldposition.s*wsize) * 20.0 + cos(worldposition.t*wsize) * 83.0)+0.1);

float wsize2 = 5.4f*1.5;
float wspeed2 = 0.2f;

float rs0a = abs(sin((worldTime*wspeed2/4.0) + (worldposition.s*wsize2) * 24.0) + 0.5);
float rs1a = abs(sin((worldTime*wspeed2/11.0) + (worldposition.t*wsize2) * 77.0  - (worldposition.z*6.0)) + 0.5);
float rs2a = abs(sin((worldTime*wspeed2/6.0) + (worldposition.s*wsize2) * 50.0 - (worldposition.t*wsize2) * 23.0) + 0.5);
float rs3a = abs(sin((worldTime*wspeed2/14.0) - (worldposition.t*wsize2) * 4.0 + (worldposition.s*wsize2) * 98.0) + 0.5);

float wsize3 = 2.0f*0.75;
float wspeed3 = 0.3f;

float rs0b = abs(sin((worldTime*wspeed3/4.0) + (worldposition.s*wsize3) * 14.0) + 0.5);
float rs1b = abs(sin((worldTime*wspeed3/11.0) + (worldposition.t*wsize3) * 37.0 + (worldposition.z*1.0)) + 0.5);
float rs2b = abs(sin((worldTime*wspeed3/6.0) + (worldposition.t*wsize3) * 47.0 - cos(worldposition.s*wsize3) * 33.0 + rs0a + rs0b) + 0.5);
float rs3b = abs(sin((worldTime*wspeed3/14.0) - (worldposition.s*wsize3) * 13.0 + sin(worldposition.t*wsize3) * 98.0 + rs0 + rs1) + 0.5);

float waves = (rs1 * rs0 + rs2 * rs3)/2.0f;
float waves2 = (rs0a * rs1a + rs2a * rs3a)/2.0f;
float waves3 = (rs0b + rs1b + rs2b + rs3b)*0.25;

float allwaves = (waves + waves2 + waves3)/3.0f;
	  allwaves *= 1.0;
	  			
				
				/*
			////////////////////////////////////RAIN WAVES////////////////////////////
			////////////////////////////////////RAIN WAVES////////////////////////////
			////////////////////////////////////RAIN WAVES////////////////////////////
float rwsize = 0.8f*3.0;
float rwspeed = 0.3f;

float r_rs0 = (sin((worldTime*rwspeed/5.0) + (worldposition.s*rwsize) * 20.0 + (worldposition.z*4.0))+0.5);
float r_rs1 = (sin((worldTime*rwspeed/7.0) + (worldposition.t*rwsize) * 27.0));
float r_rs2 = (sin((worldTime*rwspeed/2.0) + (worldposition.t*rwsize) * 60.0 - sin(worldposition.s*rwsize) * 13.0)+0.5);
float r_rs3 = (sin((worldTime*rwspeed/1.0) - (worldposition.s*rwsize) * 20.0 + cos(worldposition.t*rwsize) * 83.0)+0.5);

float rwsize2 = 0.6f*1.5;
float rwspeed2 = 0.2f;

float r_rs0a = (sin((worldTime*rwspeed2/4.0) + (worldposition.s*rwsize2) * 24.0));
float r_rs1a = (sin((worldTime*rwspeed2/11.0) + (worldposition.t*rwsize2) * 77.0  - (worldposition.z*6.0))+0.5);
float r_rs2a = (sin((worldTime*rwspeed2/6.0) + (worldposition.s*rwsize2) * 50.0 - (worldposition.t*rwsize2) * 23.0)+0.5);
float r_rs3a = (sin((worldTime*rwspeed2/14.0) - (worldposition.t*rwsize2) * 4.0 + (worldposition.s*rwsize2) * 98.0));

float rwsize3 = 0.4f*0.75;
float rwspeed3 = 0.3f;

float r_rs0b = (sin((worldTime*rwspeed3/4.0) + (worldposition.s*rwsize3) * 14.0));
float r_rs1b = (sin((worldTime*rwspeed3/11.0) + (worldposition.t*rwsize3) * 37.0 + (worldposition.z*1.0)));
float r_rs2b = (sin((worldTime*rwspeed3/6.0) + (worldposition.t*rwsize3) * 47.0 - cos(worldposition.s*rwsize3) * 33.0 + r_rs0a + r_rs0b));
float r_rs3b = (sin((worldTime*rwspeed3/14.0) - (worldposition.s*rwsize3) * 13.0 + sin(worldposition.t*rwsize3) * 98.0 + r_rs0 + r_rs1));

float rwaves = (r_rs1 * r_rs0 + r_rs2 * r_rs3)/2.0f;
float rwaves2 = (r_rs0a * r_rs1a + r_rs2a * r_rs3a)/2.0f;
float rwaves3 = (r_rs0b + r_rs1b + r_rs2b + r_rs3b)*0.25;

float rallwaves = (rwaves + rwaves2 + rwaves3)/3.0f;
	  rallwaves *= 1.0;
	  */
	  
	  
float shadingsoft = 1.0f;

	
	if (distance < drawdistance) {
		
		
		if (yDistanceSquared < drawdistancesquared) {
			

				
			if (comparedepth > 0.0f && worldposition.s < 1.0f && worldposition.s > 0.0f && worldposition.t < 1.0f && worldposition.t > 0.0f){
				//float shadowMult = min(1.0f - xzDistanceSquared / drawdistancesquared, 1.0f) * min(1.0f - yDistanceSquared / drawdistancesquared, 1.0f);
				//      shadowMult = clamp(shadowMult * 6.0f, 0.0, 1.0f);
				//shadowMult = pow(shadowMult, 0.3f);
				
				float fademult = 0.15f;
				float shadowMult = clamp((drawdistance * 0.85f * fademult) - (distance * fademult), 0.0f, 1.0f);
					 
					
					
			
					
				
				
					float zoffset = 0.00f;
					float offsetx = -0.0000f*BLURFACTOR*SHADOWOFFSET*(TimeSunset * 2.0f - 1.0f);
					float offsety = 0.0000f*BLURFACTOR*SHADOWOFFSET;
				
					
					//shadow filtering
					
					float step = 0.0f/SHADOW_RES;
					float shadowdarkness = 0.5f*SHADOW_DARKNESS;
					float diffthresh = SHADOW_CLAMP;
					float bluramount = 0.00009f*BLURFACTOR;
					
					const float confusion = 2.4f * 0.0f;

					/*
					//determine shadow depth
					float shaddepth = 0.0;
					float sds = 0.0f;
					float shaddepthspread = 1.9f * confusion;
					float stxd = -0.0010f * shaddepthspread;
					float styd = -0.0010f * shaddepthspread;
					
					for (int i = 0; i < 2; ++i) {
						stxd = -0.0010f * shaddepthspread;
						
							for (int j = 0; j < 2; ++j) {
								shaddepth =   max(shaddepth, shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx + stxd, offsety + styd) + vec2(0.0001f, 0.0001f)).z) * (256.0 - 0.05)) - zoffset, 0.0, 70.0f)/70.0f - zoffset));
								//shaddepth +=   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx + stxd, offsety + styd) + vec2(0.0001f, 0.0001f)).z) * (256.0 - 0.05)) - zoffset, 0.0, 30.0f)/30.0f - zoffset);
								stxd += 0.0005f * shaddepthspread;
								sds += 1.0f;
							}
						styd += 0.0005f * shaddepthspread;
					}
					//shaddepth /= sds;
					
					
					//fix shadow threshold
					diffthresh = 3.9f * shaddepth + 0.4f;
					
					
					//do shadows with variable blur
					shadingsharp = 1.0;
					
					int ssamp = 0;
					float shadspread = 1.9f * confusion;
					float stx = -0.0010f * shadspread;
					float sty = -0.0010f * shadspread;
					float nx = 0.0f * confusion;
					
					for (int i = 0; i < 2; ++i) {
						stx = -0.0010f * shadspread;
						
							for (int j = 0; j < 2; ++j) {
								shadingsharp +=   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx + stx * shaddepth + noiseX2*nx* shaddepth, offsety + sty * shaddepth + noiseX2*nx* shaddepth) + vec2(0.0001f, 0.0001f)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh - zoffset);
								ssamp += 1;
								stx += 0.0005f * shadspread;
							}
						sty += 0.0005f * shadspread;
					}
					
					*/
					
					
					shadingsharp =   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety) + vec2(step, step)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
					//shadingsharp =   min(shadingsharp, shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety) + vec2(-step, -step)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh - zoffset));
					
					
					
					
					//shadingsharp /= (ssamp * 1.0f);
					shadingsharp = 1.0f - shadingsharp;


					shadingsharp *= 1.0;
					shadingsharp -= 0.0;
					
					
					/*
										if (rainStrength > 0.1) {
											
											shading = 2.2;
											shadingsharp = 0.8;
										}
										*/
										
										//remove sharp shadows from water
										//shadingsharp = mix(shadingsharp, 0.2f, iswater);
										
										shading = shadingsharp;
										shading *= 0.8;
										
										//self-shadow
										shading *= texshading;
										
										
										
										//shading -= 0.2f;
										shading = clamp(shading, 0.0, 1.0);
					
					
					/////////////////////////////Skylighting///////////////////////////
					/////////////////////////////Skylighting///////////////////////////
					/////////////////////////////Skylighting///////////////////////////
					
					
									
					#ifdef SKY_LIGHTING
					
						float aospread = 5.0f;
						float trans = 0.0005 * aospread;
						float aoweight;
						float count;

						for (int i = 0; i < 5; ++i){
						
							count = i + 1;
						
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*aospread + offsetx + trans*count, noiseY2*aospread + offsety + trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*aospread + offsetx + trans*count, noiseY3*aospread + offsety - trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*aospread + offsetx - trans*count, noiseY4*aospread + offsety + trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*aospread + offsetx - trans*count, noiseY5*aospread + offsety - trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*aospread + offsetx + trans*0.0  , noiseY2*aospread + offsety + trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*aospread + offsetx + trans*0.0  , noiseY3*aospread + offsety - trans*count)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*aospread + offsetx + trans*count, noiseY4*aospread + offsety + trans*0.0  )).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							shadingao +=  shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*aospread + offsetx - trans*count, noiseY5*aospread + offsety - trans*0.0  )).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh  - zoffset)*(5.0 - i);
							
							aoweight += (5.0 - i) * 8.0f;
						
						}
						
						shadingao /= aoweight;
						
						shadingao = 1.0f - shadingao;
					
						


					
					#endif

			}
		}
	}
	
	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////
	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////
	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////

					
					vec3 npos = normalize(fragposition.xyz);

					
					vec3 specular = reflect(npos, normal);
	
					//float fresnel = distance(normal.xy, vec2(0.0f));
					//	fresnel = pow(fresnel, 6.0f);
					//	fresnel *= 3.0f;
					
					float sunlight = dot(normal, lightVector);
	

					float direct  = clamp(sin(sunlight * 3.141579f/2.0f - 0.0f) + 0.00f, 0.0, 1.00f);
						  direct = pow(direct, 0.9f);
						  //direct += max(1.0f - dot(normal, lightVector) * 3.0f - 2.0f, 0.0) * 0.05f;
						  
					float reflected = clamp(-sin(sunlight * 3.141579f/2.0f - 0.0f) + 0.95, 0.0f, 2.1f) * 0.5f;
						  reflected = pow(reflected, 3.0f);					
						  
					float ambfill = clamp(sin(sunlight * 3.141579f/2.0f - 0.0f) + 0.35, 0.0f, 1.1f);
						  ambfill = pow(ambfill, 3.0f);
						   
						  
					float spec = max(dot(specular, lightVector), 0.0f);
						  spec = pow(spec, 10.0f);
						  spec *= 2.5f;
						  spec = mix(0.0f, spec, clamp(shading, 0.0, 1.0));
						  spec = mix(0.0, spec, landx);
						  spec *= specularityDry;
						  
	
					float sunlight_direct = 0.00f + direct*1.0f;
						  //sunlight_direct = mix(1.0f, sunlight_direct, clamp(((shading-0.05)*3.0-1.1), 0.0, 1.0));
						  sunlight_direct = mix(0.0f, sunlight_direct, landx);
						  sunlight_direct = mix(sunlight_direct, 1.0f, translucent);
						  //sunlight_direct = mix(sunlight_direct, 1.0f, iswater);
						  
					float sunlight_reflected = 0.0f + reflected*1.1f;
						  //sunlight_reflected = mix(1.0f, sunlight_reflected, clamp(((shading+0.0)*1.0-0.0), 0.0, 1.0));
						  sunlight_reflected = mix(0.0f, sunlight_reflected, landx);
						  sunlight_reflected = mix(sunlight_reflected, 1.0f, translucent);
						  sunlight_reflected = mix(sunlight_reflected, 1.0f, iswater);	
						  
					float ambient_fill = 0.0f + ambfill*1.1f;
						  //ambient_fill = mix(1.0f, ambient_fill, clamp(((shading+0.0)*1.0-0.0), 0.0, 1.0));
						  ambient_fill = mix(0.0f, ambient_fill, landx);
						  ambient_fill = mix(ambient_fill, 1.0f, translucent);
						  ambient_fill = mix(ambient_fill, 1.0f, iswater);
						  
					
					shading *= sunlight_direct;
					//shading *= sunlight_reflected;
					//shading += spec;
					
					
					shading = mix(1.0, shading, landx);
				
				
				
				
				
				
				
				
				
				
				
				
				
 float gammafix = 1.0f/2.2f;
	   gammafix = mix(1.0, gammafix, landx);
 
 
//Albedo
vec4 color = texture2D(gcolor, texcoord.st);
//Linearize textures for gamma fix
	 //color *= mix(color, vec4(1.0f), 1.0f - landx);
	 color.rgb = pow(color.rgb, vec3(1.0f/gammafix));
vec3 albedo = color.rgb;

#ifdef CLAY_RENDER

color.rgb = vec3(0.5f);

#endif


	

const float rayleigh = 0.2f;

//colors for shadows/sunlight and sky
	/*
	vec3 sunrise_sun;
	 sunrise_sun.r = 1.0 * TimeSunrise;
	 sunrise_sun.g = 0.28 * TimeSunrise;
	 sunrise_sun.b = 0.0 * TimeSunrise;
	*/
	
	vec3 sunrise_sun;
	 sunrise_sun.r = 1.0 * TimeSunrise;
	 sunrise_sun.g = 0.36 * TimeSunrise;
	 sunrise_sun.b = 0.00 * TimeSunrise;
	
	vec3 sunrise_amb;
	 sunrise_amb.r = 0.00 * TimeSunrise;
	 sunrise_amb.g = 0.23 * TimeSunrise;
	 sunrise_amb.b = 0.999 * TimeSunrise;	 
	 
	
	vec3 noon_sun;
	 noon_sun.r = mix(1.00, 1.00, rayleigh) * TimeNoon;
	 noon_sun.g = mix(1.00, 0.28, rayleigh) * TimeNoon;
	 noon_sun.b = mix(0.98, 0.00, rayleigh) * TimeNoon;	 
	
	
	/*
	vec3 noon_sun;
	 noon_sun.r = 1.00 * TimeNoon;
	 noon_sun.g = 0.60 * TimeNoon;
	 noon_sun.b = 0.12 * TimeNoon;
	*/
	
	vec3 noon_amb;
	 noon_amb.r = 0.00 * TimeNoon * 1.0;
	 noon_amb.g = 0.18 * TimeNoon * 1.0;
	 noon_amb.b = 0.999 * TimeNoon * 1.0;
	
	vec3 sunset_sun;
	 sunset_sun.r = 1.0 * TimeSunset;
	 sunset_sun.g = 0.28 * TimeSunset;
	 sunset_sun.b = 0.0 * TimeSunset;
	
	vec3 sunset_amb;
	 sunset_amb.r = 0.252 * TimeSunset;
	 sunset_amb.g = 0.427 * TimeSunset;
	 sunset_amb.b = 0.999 * TimeSunset;
	
	vec3 midnight_sun;
	 midnight_sun.r = 0.3 * 0.8 * 0.05 * TimeMidnight;
	 midnight_sun.g = 0.4 * 0.8 * 0.05 * TimeMidnight;
	 midnight_sun.b = 0.8 * 0.8 * 0.05 * TimeMidnight;
	
	vec3 midnight_amb;
	 midnight_amb.r = 0.3 * 0.05 * TimeMidnight;
	 midnight_amb.g = 0.4 * 0.05 * TimeMidnight;
	 midnight_amb.b = 0.8 * 0.05 * TimeMidnight;


	vec3 sunlight_color;
	 sunlight_color.r = sunrise_sun.r + noon_sun.r + sunset_sun.r + midnight_sun.r;
	 sunlight_color.g = sunrise_sun.g + noon_sun.g + sunset_sun.g + midnight_sun.g;
	 sunlight_color.b = sunrise_sun.b + noon_sun.b + sunset_sun.b + midnight_sun.b;
	
	vec3 ambient_color;
	 ambient_color.r = sunrise_amb.r + noon_amb.r + sunset_amb.r + midnight_amb.r;
	 ambient_color.g = sunrise_amb.g + noon_amb.g + sunset_amb.g + midnight_amb.g;
	 ambient_color.b = sunrise_amb.b + noon_amb.b + sunset_amb.b + midnight_amb.b;
	 
	vec3 reflected_color;
	 reflected_color = mix(sunlight_color, ambient_color, 0.3f);
	 reflected_color = mix(vec3(0.64f, 0.73f, 0.34f), reflected_color, 0.5f);
	 reflected_color = sunlight_color;
	 
	vec3 ambfill_color;
	 ambfill_color = mix(sunlight_color, ambient_color, 0.25f);
	 
	 ambient_color = mix(ambient_color, vec3(dot(ambient_color, vec3(1.0))), SKY_DESATURATION);
	 
	 float sun_fill = 0.251f;
	
	 ambient_color = mix(ambient_color, sunlight_color, sun_fill);
	 vec3 ambient_color_rain = vec3(1.2, 1.2, 1.2) * (1.0f - TimeMidnight * 0.95f); //rain
	 ambient_color = mix(ambient_color, ambient_color_rain, rainx); //rain
	

	
		vec3 colorskyclear;
		 colorskyclear.r = ((color.r * 1.8 - 0.1) * (TimeSunrise))   +   ((color.r * 2.05 - 0.3) * (TimeNoon))   +   ((color.r * 1.8 - 0.1) * (TimeSunset))   +   (color.r * 1.0f * TimeMidnight);
		 colorskyclear.g = ((color.g * 1.8 - 0.1) * (TimeSunrise))   +   ((color.g * 2.05 - 0.4) * (TimeNoon))   +   ((color.g * 1.8 - 0.1) * (TimeSunset))   +   (color.g * 1.0f * TimeMidnight);
		 colorskyclear.b = ((color.b * 2.2 - 0.1) * (TimeSunrise))   +   ((color.b * 2.05 - 0.4) * (TimeNoon))   +   ((color.b * 2.2 - 0.1) * (TimeSunset))   +   (color.b * 1.0f * TimeMidnight);
			
			vec3 colorskyrain;
			 colorskyrain.r = ((color.r * 1.1 + 0.2) * (TimeSunrise))   +   ((color.r * 1.2 + 0.3) * (TimeNoon))   +   ((color.r * 1.1 + 0.2) * (TimeSunset))   +   (color.r * 0.1f * TimeMidnight);
			 colorskyrain.g = ((color.g * 1.1 + 0.2) * (TimeSunrise))   +   ((color.g * 1.2 + 0.3) * (TimeNoon))   +   ((color.g * 1.1 + 0.2) * (TimeSunset))   +   (color.g * 0.1f * TimeMidnight);
			 colorskyrain.b = ((color.b * 1.1 + 0.2) * (TimeSunrise))   +   ((color.b * 1.2 + 0.3) * (TimeNoon))   +   ((color.b * 1.1 + 0.2) * (TimeSunset))   +   (color.b * 0.1f * TimeMidnight);
		
			vec3 colorsky;
			 colorsky.r = mix(colorskyclear.r, colorskyrain.r, rainx);
			 colorsky.g = mix(colorskyclear.g, colorskyrain.g, rainx);
			 colorsky.b = mix(colorskyclear.b, colorskyrain.b, rainx);
			 colorsky.rgb *= 1.5f;
			
			color.r = mix(colorsky.r, color.r, landx);
			color.g = mix(colorsky.g, color.g, landx);
			color.b = mix(colorsky.b, color.b, landx);
			



//Calculate lightmap colors
sky_lightmap = pow(sky_lightmap, 3.0f);
sky_lightmap = max(sky_lightmap, 1.0 - landx);

//sky_lightmap = max(sky_lightmap, iswater);
torch_lightmap = pow(torch_lightmap, 2.5f);

 float torchwhitebalance = 0.20f;

 vec3 torchcolor;
  torchcolor.r = mix(1.00f, 1.0f, torchwhitebalance);
  torchcolor.g = mix(0.31f, 1.0f, torchwhitebalance);
  torchcolor.b = mix(0.00f, 1.0f, torchwhitebalance);
  
vec3 Specular_lightmap = vec3(spec * sunlight_color.r, spec * sunlight_color.g, spec * sunlight_color.b) * shading * (1.0f - TimeMidnight * 0.8f) * (1.0f - rainx);
	 Specular_lightmap *= pow(sky_lightmap, 0.1f);

vec3 Skylight_lightmap = vec3(sky_lightmap * ambient_color.r, sky_lightmap * ambient_color.g, sky_lightmap * ambient_color.b);

vec3 Sunlight_lightmap = vec3(shading * sunlight_color.r, shading * sunlight_color.g, shading * sunlight_color.b);
	 Sunlight_lightmap *= pow(sky_lightmap, 0.1f);
	 
vec3 Sunlight_reflected_lightmap = vec3(sunlight_reflected * reflected_color.r, sunlight_reflected * reflected_color.g, sunlight_reflected * reflected_color.b);
	 Sunlight_reflected_lightmap *= 1.5f - sky_lightmap;
	 Sunlight_reflected_lightmap *= pow(sky_lightmap, 0.2f);
	 
vec3 Sunlight_ambient_fill = vec3(ambient_fill * ambfill_color.r, ambient_fill * ambfill_color.g, ambient_fill * ambfill_color.b);
	 Sunlight_ambient_fill *= sky_lightmap;
	 
vec3 Torchlight_lightmap = vec3(torch_lightmap *  torchcolor.r, torch_lightmap *  torchcolor.g, torch_lightmap *  torchcolor.b);
	 Torchlight_lightmap.r = pow(Torchlight_lightmap.r, 1.5f);
	 Torchlight_lightmap.g = pow(Torchlight_lightmap.g, 1.5f);
	 Torchlight_lightmap.b = pow(Torchlight_lightmap.b, 1.5f);
	 
vec3 LightningFlash_lightmap = vec3(lightning_lightmap *  0.8f, lightning_lightmap *  0.7f, lightning_lightmap *  1.0f);





//RAINWET
			float dampmask = clamp(sky_lightmap * 4.0f - 1.0f, 0.0f, 1.0f) * landx * wetx;
			
			
			color.r = pow(color.r, mix(1.0f, 1.35f, dampmask));
			color.g = pow(color.g, mix(1.0f, 1.35f, dampmask));
			color.b = pow(color.b, mix(1.0f, 1.35f, dampmask));	
			
			
			
			
			
//Specular highlight


/*
	vec3 npos = normalize(fragposition.xyz);

	vec3 bump = reflect(npos, normal);
	
	float fresnel = distance(normal.xy, vec2(0.0f));
		  fresnel = pow(fresnel, 6.0f);
		  fresnel *= 3.0f;

	vec3 specularColor = vec3(sunlight_r, sunlight_g, sunlight_b) * 2.1f;
	

	float s = max(dot(normal, lightVector), 0.0);
	
	vec3 bump = specularColor * s;
		 bump *= sun_amb;
		 bump *= landx;
	*/
	
  float AO = 1.0;

#ifdef SSAO
	

  AO *= getSSAOFactor();
  
  //AO = mix(AO, 1.0f, dot(color.rgb, vec3(1.0f)) * 0.5f);
  
  AO = max(AO * 1.0f - 0.0f, 0.0f);

  //remove AO from water
  AO = mix(AO, 1.0f, iswater);
  
  //remove AO from sky
  AO = mix(1.0, AO, landx);
  
  //color.rgb *= AO;
  Sunlight_reflected_lightmap *= AO;
  Sunlight_reflected_lightmap *= AO;
  Skylight_lightmap *= AO;
  Skylight_lightmap *= AO;
  Sunlight_lightmap *= 2.0f - AO*1.0;
  Torchlight_lightmap *= AO;

#endif

float sunAOfill = 0.00f * TimeNoon + 0.000f;


//Apply different lightmaps to image
vec3 color_sky = color.rgb * (1.0f - landx);
	 color_sky = mix(color_sky, vec3(dot(color_sky, vec3(1.0))), SKY_DESATURATION);
vec3 color_skylight = color.rgb * Skylight_lightmap * landx * (1.0f - iswater) * (shadingao + SKY_LIGHTING_MIN_DARKNESS * (1.5f + (2.0f * TimeSunrise + 2.0f * TimeSunset)));
vec3 color_sunlight = color.rgb * (Sunlight_lightmap + (shadingao * sunAOfill)) * landx * (4.0f - shadingao * 3.0f) * (1.0f - iswater);
vec3 color_reflected = color.rgb * Sunlight_reflected_lightmap * landx * (shadingao)* (1.0f - iswater);
vec3 color_ambfill   = color.rgb * Sunlight_ambient_fill * landx * (shadingao)* (1.0f - iswater);

vec3 color_torchlight = color.rgb * Torchlight_lightmap * landx * (1.0f - iswater);
vec3 color_lightning = color.rgb * LightningFlash_lightmap * landx;

vec3 color_nolight = color.rgb * vec3(0.03, 0.02, 0.01);

vec3 color_water_sky = color.rgb * iswater * Skylight_lightmap;
vec3 color_water_torch = color.rgb * iswater * Torchlight_lightmap;
vec3 color_water_sunlight = color.rgb * (Sunlight_lightmap) * (iswater);

vec3 rodcolor = vec3(0.1f, 0.25f, 1.0f);


//Adjust light element levels
color_skylight         *= 1.56f * mix(1.0f, 0.15f, TimeMidnight * landx); // 0.05f
float skylight_desat    = dot(color_skylight, vec3(1.0f));
color_skylight 			= mix(color_skylight, rodcolor * skylight_desat , TimeMidnight * 0.8f * landx);

color_sunlight         *= 9.8f * mix(1.0f, 0.15f, TimeMidnight * landx);
color_sunlight         *= mix(1.0f, 0.0f, rainx); //rain
float sunlight_desat    = dot(color_sunlight, vec3(1.0f));
color_sunlight          = mix(color_sunlight, rodcolor * sunlight_desat, TimeMidnight * 0.8f);

color_reflected  *= 1.00f * mix(1.0f, 0.00f, TimeMidnight * landx) + TimeNoon * 0.333;
float color_reflected_desat = dot(color_reflected, vec3(1.0f));
color_reflected  = mix(color_reflected, color_reflected_desat * rodcolor, TimeMidnight * 0.8f * landx);
color_reflected  *= mix(1.0f, 0.0f, rainx); //rain
color_reflected  = max(color_reflected, vec3(0.0f));

color_ambfill    *= 1.090f * mix(1.0f, 0.00f, TimeMidnight * landx) + (1.0f - TimeNoon);
float color_ambfill_desat = dot(color_ambfill, vec3(1.0f));
color_ambfill 		= mix(color_ambfill, rodcolor * color_ambfill_desat, TimeMidnight * 0.8f * landx);
color_ambfill    *= mix(1.0f, 0.0f, rainx);

color_torchlight *= 7.50f;
color_lightning  *= 0.50f;

Specular_lightmap *= 05.0f;

color_sky *= 1.82f + TimeNoon * 0.45f;

color_water_torch *= 7.50f;

color_water_sky *= 1.99f;
float color_water_sky_gray = dot(color_water_sky, vec3(1.0f));
color_water_sky = mix(color_water_sky, color_water_sky_gray * rodcolor, TimeMidnight * 0.8f);

color_water_sunlight *= 0.0f;

color_nolight *= 0.09f;
float nolight_desat = dot(color_nolight, vec3(1.0f));
color_nolight = mix(color_nolight, nolight_desat * rodcolor, 0.8f);


//Add all light elements together
color.rgb = color_skylight + color_sunlight + color_reflected + color_torchlight + color_lightning + color_nolight + Specular_lightmap + color_sky + color_water_sky + color_water_torch + color_water_sunlight + color_ambfill;


//Godrays
float GRa = 0.0f;

#ifdef GODRAYS
	const float grna = 3300.0f;

	 GRa = addGodRays(0.0f, Texcoord2, noiseX3*grna, noiseX4*grna, noiseY4*grna, noiseX2*grna, noiseY2*grna)/2.0;
	 GRa = 1.0f - GRa;
	 //GRa += allwaves;
	 //GRa *= allwaves;
	 //GRa += iswater*0.25f;
#endif


color.rgb *= 0.16f;
color.b *= 1.0f;
color.rgb = mix(color.rgb * 1.3f, color.rgb * 1.0f, landx);

#ifdef PRESERVE_COLOR_RANGE
	color.rgb *= 0.5f;
#endif


//gamma fix
color.r = pow(color.r, gammafix);
color.g = pow(color.g, gammafix);
color.b = pow(color.b, gammafix);


	gl_FragData[0] = vec4(color.rgb * landx, 1.0f); //for radiosity
	//gl_FragData[0] = vec4(0.0f, 0.0f, 0.0f, 1.0f); //for optimization
	gl_FragData[0] = texture2D(gcolor, texcoord.st); //for reflections
	gl_FragData[1] = mix(texture2D(gdepth, texcoord.st), vec4(0.99999f, 0.99999f, 0.99999f, 1.0f), 1.0f - landx);
	gl_FragData[2] = vec4(normal * 0.5f + 0.5f, 1.0f);
	gl_FragData[3] = vec4(color.rgb, 1.0);
	gl_FragData[4] = texture2D(gaux1, texcoord.st);
	gl_FragData[5] = vec4(GRa, allwaves, sky_lightmap, 1.0f);
	gl_FragData[6] = texture2D(gaux3, texcoord.st);
}
