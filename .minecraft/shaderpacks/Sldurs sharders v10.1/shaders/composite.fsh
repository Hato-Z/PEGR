#version 120


//to increase shadow draw distance, edit SHADOWDISTANCE and SHADOWHPL below. Both should be equal. Needs decimal point.
//disabling is done by adding "//" to the beginning of a line.


//ADJUSTABLE VARIABLES

#define BLURFACTOR 3.5
#define SHADOW_DARKNESS 1.650   // 1.0 Is defualt darkness. 2.0 is black shadows. 0.0 is no shadows.
#define SHADOWDISTANCE 70.0 
#define SHADOW_CLAMP 0.8
#define SHADOW_RES 1024
//#define SHADOW_FILTER
	//#define VARIABLE_PENUMBRA_SHADOWS

/* SHADOWRES:1024 */
/* SHADOWHPL:70.0 */
/* gnormalFORMAT:RGB8F */

  //#define SSSM						//test
  //#define SSAO
  #define SSAO_LUMINANCE 0.0				// At what luminance will SSAO's shadows become highlights.
  #define SSAO_STRENGTH 1.60               // Too much strength causes white highlights on extruding edges and behind objects
  #define SSAO_LOOP 1						// Integer affecting samples that are taken to calculate SSAO. Higher values mean more accurate shadowing but bigger performance impact
  #define SSAO_NOISE true					// Randomize SSAO sample gathering. With noise enabled and SSAO_LOOP set to 1, you will see higher performance at the cost of fuzzy dots in shaded areas.
  #define SSAO_NOISE_AMP 0.0					// Multiplier of noise. Higher values mean SSAO takes random samples from a larger radius. Big performance hit at higher values.
  #define SSAO_MAX_DEPTH 0.9				// View distance of SSAO
  #define SSAO_SAMPLE_DELTA 0.4			// Radius of SSAO shadows. Higher values cause more performance hit.
  #define SHADOWOFFSET 0.0				// Shadow offset multiplier. Values that are too low will cause artefacts.
  //#define FXAA							// FXAA shader. Broken, but you can give it a try if you want.
  #define GODRAYS
  #define GODRAYS_EXPOSURE 0.10
  #define GODRAYS_SAMPLES 6
  #define GODRAYS_DECAY 0.95
  #define GODRAYS_DENSITY 0.65
  #define SUN_GLOW
  //#define SCREEN_SPACE_RADIOSITY			//Needed to run SSR in final.fsh

    
  //#define SKY_LIGHTING
  #define SKY_LIGHTING_SPREAD 2.63f
  #define SKY_LIGHTING_MIN_DARKNESS 0.05f 
  #define SKY_DESATURATION 0.0
  
  #define SUNLIGHT_SIZE 0.1f				//Simulates sunlight from an area instead of a point. 0.0 is a point. Higher values simulate larger sun
  #define SUNLIGHT_POWER 1.0f				//Brightness of sunlight
  #define TORCHLIGHT_POWER 0.1f
  
  //#define WATER_CAUSTICS
  
  //#define CREPUSCULAR_RAYS
  #define FOG_DENSITY 1.0
  
  //#define CLAY_RENDER

//#define PRESERVE_COLOR_RANGE

#define BANDING_FIX_FACTOR 5.5f
#define STOCHASTIC_SAMPLING

//END OF ADJUSTABLE VARIABLES






uniform sampler2D gcolor;
uniform sampler2D gdepth;
uniform sampler2D gdepthtex;
uniform sampler2D gnormal;
uniform sampler2D composite;
uniform sampler2D shadow;
uniform sampler2D watershadow;
//uniform sampler2D gaux1;
//uniform sampler2D gaux2;
//uniform sampler2D gaux3;
//uniform sampler2D gaux4;

varying vec4 texcoord;
varying vec4 lmcoord;
varying vec3 lightVector;

uniform int worldTime;

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferModelViewInverse;
uniform mat4 shadowProjection;
uniform mat4 shadowModelView;
uniform vec3 sunPosition;
uniform vec3 cameraPosition;

uniform float near;
uniform float far;
uniform float viewWidth;
uniform float viewHeight;
uniform float rainStrength;
uniform float wetness;
uniform float aspectRatio;

uniform int   isEyeInWater;
uniform float eyeAltitude;
uniform ivec2 eyeBrightness;
uniform ivec2 eyeBrightnessSmooth;



// Standard depth function.
float getDepth(vec2 coord) {
    return 2.0f * near * far / (far + near - (2.0f * texture2D(gdepth, coord).x - 1.0f) * (far - near));
}

// Convert exponential depth to linear depth
float LinearizeDepth(float depth) {
    return 2.0f * near * far / (far + near - (2.0f * depth - 1.0f) * (far - near));
}

// Convert linear depth to exponential depth
float ExponentiateDepth(float LinDepth){
	return ((( -((2.0f * near * far)/LinDepth) + far + near) / -(far - near)) + 1.0f ) / 2.0f;
}


//Auxilliary variables
//float	land 			 = texture2D(gaux1, texcoord.st).b;
//float	noblur 			 = texture2D(gaux1, texcoord.st).r;
vec3	sunPos			 = sunPosition;
vec2 	Texcoord2		 = texcoord.st;
//float 	iswater			 = texture2D(gaux1, texcoord.st).g;
vec3 	normal         	 = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
vec3 	globalNormal 	 = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
float	specularity   	 = texture2D(composite, texcoord.st).g;

//Crossfading conditionals

float rainx = clamp(rainStrength, 0.0f, 1.0f);
float wetx  = clamp(wetness, 0.0f, 1.0f);

//Lightmaps

float sky_lightmap = max(texture2D(gdepth, texcoord.st).b - 0.1f, 0.0f) * 2.0f;
float torch_lightmap = texture2D(gdepth, texcoord.st).g;
float lightning_lightmap = 0.0f;

//Standard Depth
float depth = texture2D(gdepth, texcoord.st).x;


//Calculate Time of Day

	float timefract = worldTime;
	float timePow = 3.0f;

	float TimeSunrise  = ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0) + (1.0 - (clamp(timefract, 0.0, 6000.0)/6000.0));
		  
	float TimeNoon     = ((clamp(timefract, 0.0, 6000.0)) / 6000.0) - ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0);
	  
	float TimeSunset   = ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0) - ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0);
		  
	float TimeMidnight = ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0) - ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0);


float doDistancevec2(vec2 x, vec2 y){
	return distance(x, y);
}

float doDistancevec3(vec3 x, vec3 y){
	return distance(x, y);
}
	
	
//Detect materials

	float getMaterial(vec2 coord, const int matID) {		
		float tex = texture2D(composite, coord).r;								//Call the texture carrying material
			  tex *= 255.0f; 												//Scale material info back into 0-255
			  tex = floor(tex);												//Round materials down to make sure they're integers
		
		float isMaterial;													//Create boolean for material mask
		
		if (tex == matID){
			isMaterial = 1.0f;
		} else {
			isMaterial = 0.0f;
		}
		
		return isMaterial;
	}
	
	
//Detect sky

	float getLand(vec2 coord) {		
		float tex = texture2D(composite, coord).r;								//Call the texture carrying material
			  tex *= 255.0f; 												//Scale material info back into 0-255
			  tex = floor(tex);												//Round materials down to make sure they're integers
		
		float isMaterial;													//Create boolean for material mask
		
		if (tex == 0.0f){
			isMaterial = 0.0f;
		} else {
			isMaterial = 1.0f;
		}
		
		return isMaterial;
	}


#ifdef CREPUSCULAR_RAYS
//Crepuscular Rays
	float DoAltSunGlow(float flarescale, float power){
		
			vec3 sP = sunPosition;
			
			if (sP.z > 0.0f) {
				sP.z = -sP.z;
				sP.x = -sP.x;
				sP.y = -sP.y;
				power *= 0.25f;
				}
			float c;

			vec2 lPos = sP.xy / -sP.z;
						lPos.y *= 1.39f;
						lPos.x *= 0.55f;
						
						//fix test						
						lPos.x = (lPos.x / (abs(lPos.x) + 10.0f))/(1.0/(1.0+10.0));
						lPos.y = (lPos.y / (abs(lPos.y) + 10.0f))/(1.0/(1.0+10.0));
						
						lPos = (lPos + 1.0f)/2.0f;
			
			

			vec2 flare1scale = vec2(1.7f*flarescale, 1.7f*flarescale);
			float flare1pow = 12.0f;
			vec2 flare1pos = vec2(lPos.x*aspectRatio*flare1scale.x, lPos.y*flare1scale.y);
			
			
			//float flare1 = distance(flare1pos, texcoord.st);
			float flare1 = distance(flare1pos, vec2(texcoord.s*aspectRatio*flare1scale.x, texcoord.t*flare1scale.y));
												
			
				  flare1 = 0.5 - flare1;
				  flare1 = clamp(flare1, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  //flare1 *= sunmask;
				  flare1 = pow(flare1, 10.8f);
				  
				  flare1 *= flare1pow;
				  
				  	c = flare1 * power;
		return c;			
	}

	float DrawCrepuscularRays(float noise, float landx){

		//Ray properties
			float rayCurve	  = 9.9f;
			float crepStep    = 0.45f;
			float diffthresh  = 0.1f;
			float zoffset     = 0.0f;
			float rayStep	  = 0.04f * crepStep;
			float rayStart    = 0.8f;
			
		//Initialize variables
			float crep   = 0.0f;
			float rayDepth = rayStart;
			float rayWeight = 1.0f;
			float wDepth = texture2D(gdepthtex, texcoord.st).r;
			vec2 texcoordBiased = vec2(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f);
			
			vec4 fragpositionRay = vec4(0.0f);
			vec4 worldpositionRay = vec4(0.0f);
			
		
		//Setup virtual plane on which to render shadow
			for (int i = 0; i < 100; i++){
			
				float rayDecay = i;
				
				//Abort early if ray too far
				if (rayDepth > 0.9995f) {
					break;
				}
				
				      fragpositionRay = gbufferProjectionInverse * vec4(texcoordBiased.s, texcoordBiased.t, 2.0f * min(rayDepth, 1.0f) - 1.0f, 1.0f);
					  fragpositionRay /= fragpositionRay.w;

				      worldpositionRay = gbufferModelViewInverse * fragpositionRay;
					  worldpositionRay = shadowModelView * worldpositionRay;
				
				//Shadow compare depth
				float compareRayDepth = -worldpositionRay.z;
				
					  //Shadow coordinates
					  worldpositionRay = shadowProjection * worldpositionRay;
					  worldpositionRay /= worldpositionRay.w;
					  worldpositionRay.st = worldpositionRay.st * 0.5f + 0.5f;
					  
				//Make sure ray doesn't show through geometry
				if (rayDepth < depth + (1.0f - landx)){
					
					//Draw shadow onto virtual plane
					crep += (1.0f - clamp(compareRayDepth - (0.05 + (texture2D(shadow, worldpositionRay.st).z) * (256.0 - 0.05)), 0.0, diffthresh)/(diffthresh)) * pow(rayDecay, rayCurve) * rayWeight;
										
					//Attenuate raysteps
					rayStep /= 1.0f + 0.2 * crepStep; 
					noise /= pow(rayDecay + 1, 13.0f);
					
					//Move ray forward
					rayDepth += rayStep + noise*0.4f*crepStep;
				} else {
					break;
				}
				
			}
		
		
		
		crep /= pow(100.0f, rayCurve) + 0.001f;
		
		//Anisotropy
		crep *= DoAltSunGlow(0.2f, 8000.0f) + 1.0f;
		
		//Boost in the sky
		crep *= mix(2.0f, 1.0f, landx);
		
		return crep;
		
		//return 1.0f;
	}
	
#endif

// Alternate projected depth (used by SSAO, probably AA too)
float getProDepth(vec2 coord) {
	float depth = texture2D(gdepth, coord).x;
	return ( 2.0f * near ) / ( far + near - depth * ( far - near ) );
}
	
	
//Noise
float noise(vec2 coord, float width, float height) { //generating noise/pattern texture for dithering
  float noiseX = ((fract(1.0f-coord.s*(width/2.0f))*0.25f)+(fract(coord.t*(height/2.0f))*0.75f))*2.0f-1.0f;
  float noiseY = ((fract(1.0f-coord.s*(width/2.0f))*0.75f)+(fract(coord.t*(height/2.0f))*0.25f))*2.0f-1.0f;

  //generate SSAO noise
  noiseX = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f))) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  noiseY = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f)*2.0f)) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  
  return noiseX;
}

#ifdef SSAO

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
	
	
	ao /= 16.0f;
	
	ao = 1.0f-ao;	
  ao = clamp(ao, 0.0f, 0.5f) * 2.0f;
	
  return ao;
}

#endif






#ifdef SSSM

float diffareaS = 0.1f; //self-shadowing reduction
float gdisplaceS = 0.30f; //gauss bell center


float compareDepthsS(in float depth1, in float depth2, in int zfar) {  
  float garea = 0.1f; //gauss bell width    
  float diff = (depth1 - depth2) * 100.0f; //depth difference (0-100)
  //reduce left bell width to avoid self-shadowing 
  
  if (diff < gdisplaceS) {
    garea = diffareaS;
  } else {
    zfar = 1;
  }


  float gauss = pow(2.7182f,-2.0f*(diff-gdisplaceS)*(diff-gdisplaceS)/(garea*garea));
  return gauss;
} 



float calSSSM(float depth, float dw, float dh) {  
  float temp = 0.0f;
  float temp2 = 0.0f;
  dw *= 2.0f;
  dh *= 2.0f;
  float coordw = texcoord.x + dw/(depth*0.2f + 0.1f);
  float coordh = texcoord.y + dh/(depth*0.2f + 0.1f);
  float coordw2 = texcoord.x - dw/(depth*0.2f + 0.1f);
  float coordh2 = texcoord.y - dh/(depth*0.2f + 0.1f);

  if (coordw  < 1.0f && coordw  > 0.0f && coordh < 1.0f && coordh  > 0.0f){
    vec2 coord = vec2(coordw , coordh);
    vec2 coord2 = vec2(coordw2, coordh2);
    int zfar = 0;
    temp = compareDepthsS(depth, getProDepth(coord),zfar);
  }

  return temp;  
}  



float getSSSM() {

	float depth = getProDepth(texcoord.st);
	float cdepth = texture2D(gdepth, texcoord.st).x;
	
	float shad = 0.0f;
	
	vec2 shadvector = lightVector.xy * 0.00025f;
		 shadvector = shadvector * vec2(1.0f, aspectRatio);
	
		for (int i = 1; i < 12; ++i) {
			shad = max(shad, calSSSM(depth, shadvector.x * i, shadvector.y * i));
		}
		
		shad /= 1.0f;
		
		shad = 1.0f - shad;
		
	return shad;
  
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
			//lightPos.y *= 1.39f;
			//lightPos.x *= 0.76f;
			//lightPos.y *= aspectRatio;
			lightPos.x *= 1.0f/aspectRatio;
			lightPos.xy *= 1.40f;
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

					sample = 1.0f - getLand(tx + vec2(noise*delta.x, noise*delta.y));
					sample += 1.0f - getLand(tx + vec2(noise2*delta.x, noise2*delta.y));
					sample += 1.0f - getLand(tx + vec2(noise3*delta.x, noise3*delta.y));
					sample += 1.0f - getLand(tx + vec2(noise4*delta.x, noise4*delta.y));
					sample += 1.0f - getLand(tx + vec2(noise5*delta.x, noise5*delta.y));
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



#ifdef SUN_GLOW
	
	vec3 DoSunGlow(float flarescale, vec3 suncol, float power, float landx){
			vec3 sP = sunPosition;
			vec3 c;

			vec2 lPos = sP.xy / -sP.z;
						lPos.y *= 1.39f;
						lPos.x *= 0.55f;
						
						//fix test						
						lPos.x = (lPos.x / (abs(lPos.x) + 10.0f))/(1.0/(1.0+10.0));
						lPos.y = (lPos.y / (abs(lPos.y) + 10.0f))/(1.0/(1.0+10.0));
						
						lPos = (lPos + 1.0f)/2.0f;
			
			

			vec2 flare1scale = vec2(1.7f*flarescale, 1.7f*flarescale);
			float flare1pow = 12.0f;
			vec2 flare1pos = vec2(lPos.x*aspectRatio*flare1scale.x, lPos.y*flare1scale.y);
			
			
			//float flare1 = distance(flare1pos, texcoord.st);
			float flare1 = distance(flare1pos, vec2(texcoord.s*aspectRatio*flare1scale.x, texcoord.t*flare1scale.y));
												
			
				  flare1 = 0.5 - flare1;
				  flare1 = clamp(flare1, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  //flare1 *= sunmask;
				  flare1 = pow(flare1, 5.8f);
				  
				  flare1 *= flare1pow;
				  
				  	c = flare1 * suncol * power * (1.0f - landx);
		return c;			
	}

#endif





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void main() {

/////////////////////////GET MATERIALS/////////////////////
/////////////////////////GET MATERIALS/////////////////////
/////////////////////////GET MATERIALS/////////////////////
float land 			= getLand(texcoord.st);
float landx 		= land;
float grass			= getMaterial(texcoord.st, 2);
float leaves 		= getMaterial(texcoord.st, 3);
float ice			= getMaterial(texcoord.st, 4);
float hand			= getMaterial(texcoord.st, 5);
float translucent   = getMaterial(texcoord.st, 6);
float iswater = 0.0f;
if (texture2D(composite, texcoord.st).b >= 0.25f) {
	iswater = 1.0f;
}



//Curve times
//Curve times
//Curve times
		  TimeSunrise  = pow(TimeSunrise, timePow);
		  TimeNoon     = pow(TimeNoon, 1.0f/timePow);
		  TimeSunset   = pow(TimeSunset, timePow);
		  TimeMidnight = pow(TimeMidnight, 1.0f/timePow);

float noiseamp = 0.3f;
				
				#ifdef STOCHASTIC_SAMPLING
					
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
						
					#ifdef SKY_LIGHTING
						float width5 = 4.0f;
						float height5 = 4.0f;
						float noiseX5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.25f)+(fract(Texcoord2.t*(height5/2.0f))*0.75f))*2.0f-1.0f;
						float noiseY5 = ((fract(1.0f-Texcoord2.s*(width5/2.0f))*0.75f)+(fract(Texcoord2.t*(height5/2.0f))*0.25f))*2.0f-1.0f;

						
							noiseX5 = clamp(fract(sin(dot(Texcoord2 ,vec2(11.9898f,68.633f))) * 21178.5453f),0.0f,1.0f)*2.0f-1.0f;
							noiseY5 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898f,71.233f)*2.0f)) * 6958.5453f),0.0f,1.0f)*2.0f-1.0f;
						
						noiseX5 *= (0.0005f*noiseamp);
						noiseY5 *= (0.0005f*noiseamp);
					#endif	
					
				#else
					float noiseX2 = 0.0f;
					float noiseY2 = 0.0f;
					float noiseX3 = 0.0f;
					float noiseY3 = 0.0f;
					float noiseX4 = 0.0f;
					float noiseY4 = 0.0f;
					float noiseX5 = 0.0f;
					float noiseY5 = 0.0f;
				#endif
//

	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * texture2D(gdepth, texcoord.st).x - 1.0f, 1.0f);
	vec4 fragpositionw = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * texture2D(gdepthtex, texcoord.st).x - 1.0f, 1.0f);
	fragposition /= fragposition.w;
	fragpositionw /= fragpositionw.w;
	
	#ifdef SHADOWDISTANCE
	float drawdistance = SHADOWDISTANCE;
	float drawdistancesquared = pow(drawdistance, 2.0f);
	#endif
	
	float distance = sqrt(fragposition.x * fragposition.x + fragposition.y * fragposition.y + fragposition.z * fragposition.z);

	float shading = 1.0f;
	float shadingsharp = 1.0f;
	float shadingao = 1.0f;
	
	
	vec4 worldposition = vec4(0.0);
	vec4 worldpositionw = vec4(0.0);
			
	worldposition = gbufferModelViewInverse * fragposition;	
	worldpositionw = gbufferModelViewInverse * fragpositionw;	
	
	float xzDistanceSquared = worldposition.x * worldposition.x + worldposition.z * worldposition.z;
	float yDistanceSquared  = worldposition.y * worldposition.y;

	
			worldposition = shadowModelView * worldposition;
			worldpositionw = shadowModelView * worldpositionw;
			float comparedepth = -worldposition.z;
			float comparedepthWater = -worldpositionw.z;
			worldposition = shadowProjection * worldposition;
	vec4 worldpositionWaves = worldpositionw;
			worldpositionw = shadowProjection * worldpositionw;
			worldposition /= worldposition.w;
			worldpositionw /= worldpositionw.w;
			
			worldposition.st = worldposition.st * 0.5f + 0.5f;
			worldpositionw.st = worldpositionw.st * 0.5f + 0.5f;
			
			
			
			
			////////////////////////////////////WAVES////////////////////////////
			////////////////////////////////////WAVES////////////////////////////
			////////////////////////////////////WAVES////////////////////////////
float wsize = 0.1f*1.5;
float wspeed = 0.3f;

float rs0 = abs(sin((worldTime*wspeed/5.0) + (worldpositionWaves.s*wsize) * 20.0 + (worldpositionWaves.z*4.0))+0.2);
float rs1 = abs(sin((worldTime*wspeed/7.0) + (worldpositionWaves.t*wsize) * 27.0) + 0.5);
float rs2 = abs(sin((worldTime*wspeed/2.0) + (worldpositionWaves.t*wsize) * 60.0 - sin(worldpositionWaves.s*wsize) * 13.0)+0.4);
float rs3 = abs(sin((worldTime*wspeed/1.0) - (worldpositionWaves.s*wsize) * 20.0 + cos(worldpositionWaves.t*wsize) * 83.0)+0.1);

float wsize2 = 0.05f*0.75;
float wspeed2 = 0.2f;

float rs0a = abs(sin((worldTime*wspeed2/4.0) + (worldpositionWaves.s*wsize2) * 24.0) + 0.5);
float rs1a = abs(sin((worldTime*wspeed2/11.0) + (worldpositionWaves.t*wsize2) * 77.0  - (worldpositionWaves.z*6.0)) + 0.5);
float rs2a = abs(sin((worldTime*wspeed2/6.0) + (worldpositionWaves.s*wsize2) * 50.0 - (worldpositionWaves.t*wsize2) * 23.0) + 0.5);
float rs3a = abs(sin((worldTime*wspeed2/14.0) - (worldpositionWaves.t*wsize2) * 4.0 + (worldpositionWaves.s*wsize2) * 98.0) + 0.5);

float wsize3 = 0.01f*0.125;
float wspeed3 = 0.3f;

float rs0b = abs(sin((worldTime*wspeed3/4.0) + (worldpositionWaves.s*wsize3) * 14.0) + 0.5);
float rs1b = abs(sin((worldTime*wspeed3/11.0) + (worldpositionWaves.t*wsize3) * 37.0 + (worldpositionWaves.z*1.0)) + 0.5);
float rs2b = abs(sin((worldTime*wspeed3/6.0) + (worldpositionWaves.t*wsize3) * 47.0 - cos(worldpositionWaves.s*wsize3) * 33.0 + rs0a + rs0b) + 0.5);
float rs3b = abs(sin((worldTime*wspeed3/14.0) - (worldpositionWaves.s*wsize3) * 13.0 + sin(worldpositionWaves.t*wsize3) * 98.0 + rs0 + rs1) + 0.5);

float waves = (rs1 * rs0 + rs2 * rs3)/2.0f;
float waves2 = (rs0a * rs1a + rs2a * rs3a)/2.0f;
float waves3 = (rs0b + rs1b + rs2b + rs3b)*0.25;

float allwaves = (waves + waves2 + waves3)/3.0f;
	  allwaves *= 1.0;			
			
			
			
			////////////////////////////////////CAUSTIC WAVES////////////////////////////
			////////////////////////////////////CAUSTIC WAVES////////////////////////////
			////////////////////////////////////CAUSTIC WAVES////////////////////////////
float wsizeC = 9.0f*3.0;
float wspeedC = 0.3f;

float rs0C = abs(sin((worldTime*wspeed/5.0) + (worldposition.s*wsizeC) * 20.0 + (worldposition.z*4.0))+0.2);
float rs1C = abs(sin((worldTime*wspeedC/7.0) + (worldposition.t*wsizeC) * 27.0) + 0.5);
float rs2C = abs(sin((worldTime*wspeedC/2.0) + (worldposition.t*wsizeC) * 60.0 - sin(worldposition.s*wsizeC) * 13.0)+0.4);
float rs3C = abs(sin((worldTime*wspeedC/1.0) - (worldposition.s*wsizeC) * 20.0 + cos(worldposition.t*wsizeC) * 83.0)+0.1);

float wsizeC2 = 5.4f*1.5;
float wspeedC2 = 0.2f;

float rs0Ca = abs(sin((worldTime*wspeedC2/4.0) + (worldposition.s*wsizeC2) * 24.0) + 0.5);
float rs1Ca = abs(sin((worldTime*wspeedC2/11.0) + (worldposition.t*wsizeC2) * 77.0  - (worldposition.z*6.0)) + 0.5);
float rs2Ca = abs(sin((worldTime*wspeedC2/6.0) + (worldposition.s*wsizeC2) * 50.0 - (worldposition.t*wsizeC2) * 23.0) + 0.5);
float rs3Ca = abs(sin((worldTime*wspeedC2/14.0) - (worldposition.t*wsizeC2) * 4.0 + (worldposition.s*wsizeC2) * 98.0) + 0.5);

float wsizeC3 = 2.0f*0.75;
float wspeedC3 = 0.3f;

float rs0Cb = abs(sin((worldTime*wspeedC3/4.0) + (worldposition.s*wsizeC3) * 14.0) + 0.5);
float rs1Cb = abs(sin((worldTime*wspeedC3/11.0) + (worldposition.t*wsizeC3) * 37.0 + (worldposition.z*1.0)) + 0.5);
float rs2Cb = abs(sin((worldTime*wspeedC3/6.0) + (worldposition.t*wsizeC3) * 47.0 - cos(worldposition.s*wsizeC3) * 33.0 + rs0Ca + rs0Cb) + 0.5);
float rs3Cb = abs(sin((worldTime*wspeedC3/14.0) - (worldposition.s*wsizeC3) * 13.0 + sin(worldposition.t*wsizeC3) * 98.0 + rs0C + rs1C) + 0.5);

float wavesC = (rs1C * rs0C + rs2C * rs3C)/2.0f;
float wavesC2 = (rs0Ca * rs1Ca + rs2Ca * rs3Ca)/2.0f;
float wavesC3 = (rs0Cb + rs1Cb + rs2Cb + rs3Cb)*0.25;

float allwavesC = (wavesC + wavesC2 + wavesC3)/3.0f;
	  allwavesC *= 1.0;
	  allwavesC = wavesC + wavesC2;
	  allwavesC /= 2.0f;
	  			
				
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

float shadowMult = 0.0f;
	  
float shadingsoft = 1.0f;

float shadingWater = 1.0f;

	
	if (distance < drawdistance) {
		
		
		if (yDistanceSquared < drawdistancesquared) {
			

				
			if (comparedepth > 0.0f && worldposition.s < 1.0f && worldposition.s > 0.0f && worldposition.t < 1.0f && worldposition.t > 0.0f){
				//float shadowMult = min(1.0f - xzDistanceSquared / drawdistancesquared, 1.0f) * min(1.0f - yDistanceSquared / drawdistancesquared, 1.0f);
				//      shadowMult = clamp(shadowMult * 6.0f, 0.0, 1.0f);
				//shadowMult = pow(shadowMult, 0.3f);
				
				float fademult = 0.15f;
					  shadowMult = clamp((drawdistance * 0.85f * fademult) - (distance * fademult), 0.0f, 1.0f);
					 
					
					
			
					
				
				
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
					
					
					#ifdef SHADOW_FILTER
						
						float vpsdepth = 0.0f;
						float vpsconfusion = 3.5f * SUNLIGHT_SIZE;
					
						#ifdef VARIABLE_PENUMBRA_SHADOWS
							
							float maxCompareDepth = 60.0f;
							vec2 vpsSpread = vec2(5.8f/2048.0f) * vpsconfusion;
							vec2 vpsc = vec2(-2.0f);
							float vpssamp = 0.0f;
							float vpscurve = 0.45f;
							
							for(int i = 0; i < 5; ++i){
								vpsc.x = -2.0f;
									
								for(int j = 0; j < 5; ++j){
									vpsdepth += pow((clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety) + vpsSpread * vpsc).z) * (256.0 - 0.05)), 0.0, maxCompareDepth)/(maxCompareDepth)), 1.0f/vpscurve);
									vpsc.x += 1.0f;
									vpssamp += 1.0f;
								}
								vpsc.y += 1.0f;
							}
							
							vpsdepth = pow(vpsdepth, vpscurve);
							vpsdepth /= vpssamp;
							vpsdepth *= vpsconfusion;
							
						#endif
						
						diffthresh = diffthresh * vpsdepth * 45.0f + diffthresh;
					
						float sfx = 20.8f * vpsdepth/2048.0f;
						float sfy = 20.8f * vpsdepth/2048.0f;
						
						float sfxc = -5.0f;
						float sfyc = -5.0f;
						float sfsamp = 0.0f;
					
						for(int i = 0; i < 11; ++i){
						
							//break early if not sufficient blur
							if(vpsdepth < 0.001f){
								shadingsharp = shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
								sfsamp = 1.0f;
								break;
							}
						
							sfxc = -5.0f;
						
							for(int j = 0; j < 11; ++j){
							
								shadingsharp += shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety) + vec2(sfx * sfxc, sfy * sfyc)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
								sfxc += 1.0f;
								sfsamp += 1.0f;
							
							}
							
							sfyc += 1.0f;
							
						}
						
						shadingsharp /= sfsamp;
					
					
					#else
	
					shadingsharp =   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety) + vec2(step, step)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
					
					#endif
					
					//Water sunlight occlusion
					shadingWater =   (clamp(comparedepthWater - (0.05 + (texture2D(shadow, worldpositionw.st + vec2(offsetx, offsety) + vec2(step, step)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
					shadingWater *= iswater;
					shadingWater = 1.0f - shadingWater;
					
					
					
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
										//shading *= texshading;
										
										
										
										//shading -= 0.2f;
										shading = clamp(shading, 0.0, 1.0);
					
					
					#ifdef WATER_CAUSTICS
						float wshadow = 1.0f;
						
						float wdepth = 1.0f;
						float wdepththresh = 20.0f;
							wdepth =   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(watershadow, worldposition.st + vec2(offsetx, offsety) + vec2(step, step)).z) * (256.0 - 0.05)) - zoffset, 0.0, wdepththresh)/(wdepththresh) - zoffset);
							wdepth = pow(wdepth, 0.25f);

						//caustics
							float caustics = pow(abs(sin(abs(allwavesC)*3.1415f)), 10.0f * max(1.0f - wdepth * 1.2f, 0.0f) + 1.0f) * 5.9f;
								  caustics = mix(caustics, 1.0f, pow(wdepth, 1.9f))*wdepth*1.0f;
								  caustics *= sky_lightmap;
							
							float wdiffthresh = 0.9f;

							
							wshadow =   shadowMult * (clamp(comparedepth - (0.05 + (texture2D(watershadow, worldposition.st + vec2(offsetx, offsety) + vec2(step, step)).z) * (256.0 - 0.05)) - zoffset, 0.0, wdiffthresh)/(wdiffthresh) - zoffset);
							wshadow = 1.0f - wshadow;
							
							wshadow = mix(wshadow + caustics * (max(1.0f - wdepth, 0.0f)) * sky_lightmap, 1.0f, wshadow);
							shading *= wshadow;
							
					#endif
					
					
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

					
					vec3 halfVector = normalize(lightVector - npos);
					float specular 	= max(0.0f, dot(halfVector, normal));
	
					
					float sunlight = dot(normal, lightVector);
					float sunlightback = dot(normal, -lightVector);
	

					float direct  = clamp(sin(sunlight * 3.141579f/2.0f - 0.0f) + SUNLIGHT_SIZE/2.0f, 0.0, 1.00f + SUNLIGHT_SIZE/2.0f) / (1.0f + SUNLIGHT_SIZE/2.0f);
						  direct = pow(direct, 1.0f);
						  //direct += max(1.0f - dot(normal, lightVector) * 3.0f - 2.0f, 0.0) * 0.05f;
					
					float directleaves = clamp(sin(sunlight * 3.141579f/2.0f - 0.0f) + SUNLIGHT_SIZE/2.0f, 0.0, 1.00f + SUNLIGHT_SIZE/2.0f) / (1.0f + SUNLIGHT_SIZE/2.0f);
						  directleaves += (clamp(sin(sunlightback * 3.141579f/2.0f - 0.0f) + SUNLIGHT_SIZE/2.0f, 0.0, 1.00f + SUNLIGHT_SIZE/2.0f) / (1.0f + SUNLIGHT_SIZE/2.0f)) * 0.75f;
						 
						  //directleaves = 1.0f;
						  //direct += max(1.0f - dot(normal, lightVector) * 3.0f - 2.0f, 0.0) * 0.05f;
						  
					float reflected = clamp(-sin(sunlight * 3.141579f/2.0f - 0.0f) + 0.95, 0.0f, 2.1f) * 0.5f;
						  reflected = pow(reflected, 3.0f);					
						  
					float ambfill = clamp(sin(sunlight * 3.141579f/2.0f - 0.0f) + 0.55, 0.0f, 1.5f);
						  ambfill = pow(ambfill, 3.0f);
						   
						  
					float spec = specular;
						  spec = pow(spec, 50.0f);
						  spec *= 5.0f;
						  spec = mix(0.0f, spec, clamp(shading, 0.0, 1.0));
						  spec = mix(0.0, spec, landx);
						  spec *= specularity;
						  spec = mix(0.0f, spec, shadowMult);
						  
	
					float sunlight_direct = 0.00f + mix(direct, directleaves, leaves);
						  //sunlight_direct = mix(1.0f, sunlight_direct, clamp(((shading-0.05)*3.0-1.1), 0.0, 1.0));
						  sunlight_direct = mix(0.0f, sunlight_direct, landx);
						  sunlight_direct = mix(sunlight_direct, 1.0f, grass);
						  sunlight_direct = mix(sunlight_direct, 1.0f, iswater);
						  //sunlight_direct = mix(sunlight_direct, 1.0f, iswater);
						  //sunlight_direct = mix(1.0f, sunlight_direct, shadowMult);
						  
						  
					float sunlight_reflected = 0.0f + reflected*1.1f;
						  //sunlight_reflected = mix(1.0f, sunlight_reflected, clamp(((shading+0.0)*1.0-0.0), 0.0, 1.0));
						  sunlight_reflected = mix(0.0f, sunlight_reflected, landx);
						  sunlight_reflected = mix(sunlight_reflected, 1.0f, grass);
						  //sunlight_reflected = mix(sunlight_reflected, 1.0f, iswater);
						  sunlight_reflected = mix(0.0f, sunlight_reflected, shadowMult);
						  
					float ambient_fill = 0.0f + ambfill*1.1f;
						  //ambient_fill = mix(1.0f, ambient_fill, clamp(((shading+0.0)*1.0-0.0), 0.0, 1.0));
						  ambient_fill = mix(0.0f, ambient_fill, landx);
						  ambient_fill = mix(ambient_fill, 1.0f, grass);
						  //ambient_fill = mix(ambient_fill, 1.0f, iswater);
						  ambient_fill = mix(0.0f, ambient_fill, shadowMult);
						  
						  
						  

					
					
					
					
					
					
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
	 
	 //Additional texture fix
	 //color.rgb = pow(color.rgb, vec3(1.25f));
vec3 albedo = color.rgb;

#ifdef CLAY_RENDER

color.rgb = vec3(0.5f);

#endif


	
const float rayleigh = 0.0f;


//colors for shadows/sunlight and sky
	/*
	vec3 sunrise_sun;
	 sunrise_sun.r = 1.0 * TimeSunrise;
	 sunrise_sun.g = 0.28 * TimeSunrise;
	 sunrise_sun.b = 0.0 * TimeSunrise;
	*/
	
	vec3 sunrise_sun;
	 sunrise_sun.r = 1.0 * TimeSunrise;
	 sunrise_sun.g = 0.46 * TimeSunrise;
	 sunrise_sun.b = 0.00 * TimeSunrise;
	 sunrise_sun *= 0.85f;
	
	vec3 sunrise_amb;
	 sunrise_amb.r = 0.00 * TimeSunrise;
	 sunrise_amb.g = 0.23 * TimeSunrise;
	 sunrise_amb.b = 0.999 * TimeSunrise;	
	 sunrise_amb = mix(sunrise_amb, vec3(1.0f), 0.2f);
	 
	
	vec3 noon_sun;
	 noon_sun.r = mix(1.00, 1.00, rayleigh) * TimeNoon;
	 noon_sun.g = mix(1.00, 0.48, rayleigh) * TimeNoon;
	 noon_sun.b = mix(1.00, 0.00, rayleigh) * TimeNoon;	 
	
	
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
	 sunset_sun.g = 0.38 * TimeSunset;
	 sunset_sun.b = 0.0 * TimeSunset;
	 sunset_sun *= 0.85f;
	
	vec3 sunset_amb;
	 sunset_amb.r = 0.252 * TimeSunset;
	 sunset_amb.g = 0.427 * TimeSunset;
	 sunset_amb.b = 0.999 * TimeSunset;
	
	vec3 midnight_sun;
	 midnight_sun.r = 0.45 * 0.8 * 0.25 * TimeMidnight;
	 midnight_sun.g = 0.6 * 0.8 * 0.25 * TimeMidnight;
	 midnight_sun.b = 0.8 * 0.8 * 0.25 * TimeMidnight;
	
	vec3 midnight_amb;
	 midnight_amb.r = 0.3 * 0.15 * TimeMidnight;
	 midnight_amb.g = 0.4 * 0.15 * TimeMidnight;
	 midnight_amb.b = 0.8 * 0.15 * TimeMidnight;


	vec3 sunlight_color;
	 sunlight_color.r = sunrise_sun.r + noon_sun.r + sunset_sun.r + midnight_sun.r;
	 sunlight_color.g = sunrise_sun.g + noon_sun.g + sunset_sun.g + midnight_sun.g;
	 sunlight_color.b = sunrise_sun.b + noon_sun.b + sunset_sun.b + midnight_sun.b;
	
	vec3 ambient_color;
	 ambient_color.r = sunrise_amb.r + noon_amb.r + sunset_amb.r + midnight_amb.r;
	 ambient_color.g = sunrise_amb.g + noon_amb.g + sunset_amb.g + midnight_amb.g;
	 ambient_color.b = sunrise_amb.b + noon_amb.b + sunset_amb.b + midnight_amb.b;
	 
	vec3 reflected_color;
	 reflected_color = mix(sunlight_color, ambient_color, 0.25f);
	 //reflected_color = mix(vec3(0.64f, 0.73f, 0.34f), reflected_color, 0.5f);
	 //reflected_color = sunlight_color;
	 
	vec3 ambfill_color;
	 ambfill_color = mix(sunlight_color, ambient_color, 0.55f);
	 
	 ambient_color = mix(ambient_color, vec3(dot(ambient_color, vec3(1.0))), SKY_DESATURATION);
	 
	vec3 skycolor = mix(sunlight_color, vec3(1.0f), 0.99f);
	 
	 float sun_fill = 0.201f;
	
	 ambient_color = mix(ambient_color, sunlight_color, sun_fill);
	 vec3 ambient_color_rain = vec3(1.78, 1.78, 1.78) * (1.0f - TimeMidnight * 0.95f); //rain
	 ambient_color = mix(ambient_color, ambient_color_rain, rainx); //rain
	

	
		vec3 colorskyclear;
		 colorskyclear.r = ((color.r * 1.8 - 0.1) * (TimeSunrise))   +   ((color.r * 2.05 - 0.4) * (TimeNoon))   +   ((color.r * 1.8 - 0.1) * (TimeSunset))   +   (color.r * 1.0f * TimeMidnight);
		 colorskyclear.g = ((color.g * 1.8 - 0.1) * (TimeSunrise))   +   ((color.g * 2.05 - 0.4) * (TimeNoon))   +   ((color.g * 1.8 - 0.1) * (TimeSunset))   +   (color.g * 1.0f * TimeMidnight);
		 colorskyclear.b = ((color.b * 2.2 - 0.1) * (TimeSunrise))   +   ((color.b * 2.05 - 0.4) * (TimeNoon))   +   ((color.b * 2.2 - 0.1) * (TimeSunset))   +   (color.b * 1.0f * TimeMidnight);
			
			vec3 colorskyrain;
			 colorskyrain.r = ((color.r * 1.4 + 0.2) * (TimeSunrise))   +   ((color.r * 1.4 + 0.3) * (TimeNoon))   +   ((color.r * 1.3 + 0.2) * (TimeSunset))   +   (color.r * 0.1f * TimeMidnight);
			 colorskyrain.g = ((color.g * 1.4 + 0.2) * (TimeSunrise))   +   ((color.g * 1.4 + 0.3) * (TimeNoon))   +   ((color.g * 1.3 + 0.2) * (TimeSunset))   +   (color.g * 0.1f * TimeMidnight);
			 colorskyrain.b = ((color.b * 1.4 + 0.2) * (TimeSunrise))   +   ((color.b * 1.4 + 0.3) * (TimeNoon))   +   ((color.b * 1.3 + 0.2) * (TimeSunset))   +   (color.b * 0.1f * TimeMidnight);
		
			vec3 colorsky;
			 colorsky.r = mix(colorskyclear.r, colorskyrain.r, rainx);
			 colorsky.g = mix(colorskyclear.g, colorskyrain.g, rainx);
			 colorsky.b = mix(colorskyclear.b, colorskyrain.b, rainx);
			 colorsky.rgb *= 1.5f;
			
			color.r = mix(colorsky.r, color.r, landx);
			color.g = mix(colorsky.g, color.g, landx);
			color.b = mix(colorsky.b, color.b, landx);
			
	//Saturate sunlight colors
	sunlight_color = pow(sunlight_color, vec3(3.0f));
			



//Calculate lightmap colors
sky_lightmap = pow(sky_lightmap, 2.9f);
sky_lightmap = max(sky_lightmap, 1.0 - landx);

//sky_lightmap = max(sky_lightmap, iswater);
torch_lightmap = pow(torch_lightmap, 3.0f);

 float torchwhitebalance = 0.10f;

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
	 Sunlight_reflected_lightmap *= pow(sky_lightmap, 0.8f);
	 
vec3 Sunlight_ambient_fill = vec3(ambient_fill * ambfill_color.r, ambient_fill * ambfill_color.g, ambient_fill * ambfill_color.b);
	 Sunlight_ambient_fill *= sky_lightmap;
	 
vec3 Torchlight_lightmap = vec3(torch_lightmap *  torchcolor.r, torch_lightmap *  torchcolor.g, torch_lightmap *  torchcolor.b);
	 Torchlight_lightmap.r = pow(Torchlight_lightmap.r, 1.1f);
	 Torchlight_lightmap.g = pow(Torchlight_lightmap.g, 1.1f);
	 Torchlight_lightmap.b = pow(Torchlight_lightmap.b, 1.1f);
	 
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
	//AO = mix(AO, 1.0f, iswater);
  
  //remove AO from sky
  AO = mix(1.0, AO, landx);
  
  //color.rgb *= AO;
  Sunlight_reflected_lightmap *= AO;
  Sunlight_reflected_lightmap *= AO;
  //Sunlight_reflected_lightmap *= AO;
  Sunlight_ambient_fill *= AO;
  Sunlight_ambient_fill *= AO;
  Sunlight_ambient_fill *= AO;
  Skylight_lightmap *= AO;
  //Skylight_lightmap *= AO;
  Sunlight_lightmap *= 2.5f - AO*1.5;
  //Sunlight_lightmap *= AO;
  Torchlight_lightmap *= AO;

#endif




float screenshad = 1.0f;

#ifdef SSSM

	screenshad = getSSSM();
	
	Sunlight_lightmap *= screenshad;


//test


#endif

float sunAOfill = 0.00f * TimeNoon + 0.000f;


//Apply different lightmaps to image
vec3 color_sky = color.rgb * (1.0f - landx);
	 color_sky = mix(color_sky, vec3(dot(color_sky, vec3(1.0))), SKY_DESATURATION);
	 color_sky *= skycolor;
vec3 color_skylight = color.rgb * Skylight_lightmap * landx * (shadingao + SKY_LIGHTING_MIN_DARKNESS * (1.5f + (2.0f * TimeSunrise + 2.0f * TimeSunset)));
vec3 color_sunlight = color.rgb * (Sunlight_lightmap + (shadingao * sunAOfill)) * landx * (4.0f - shadingao * 3.0f);
vec3 color_reflected = color.rgb * Sunlight_reflected_lightmap * landx * (shadingao);
vec3 color_ambfill   = color.rgb * Sunlight_ambient_fill * landx * (shadingao);

vec3 color_torchlight = color.rgb * Torchlight_lightmap * landx;
vec3 color_lightning = color.rgb * LightningFlash_lightmap * landx;

vec3 color_nolight = color.rgb * vec3(0.03, 0.02, 0.01);

vec3 color_water_sky = color.rgb * iswater * Skylight_lightmap;
vec3 color_water_torch = color.rgb * iswater * Torchlight_lightmap;
vec3 color_water_sunlight = color.rgb * (Sunlight_lightmap) * (iswater);

vec3 rodcolor = vec3(0.1f, 0.25f, 1.0f);


//Adjust light element levels
#ifdef SKY_LIGHTING
color_skylight         *= 4.96f * mix(1.0f, 0.15f, TimeMidnight * landx); // 0.05f
#else
color_skylight         *= 3.36f * mix(1.0f, 0.15f, TimeMidnight * landx); // 0.05f
#endif
float skylight_desat    = dot(color_skylight, vec3(1.0f));
color_skylight 			= mix(color_skylight, rodcolor * skylight_desat , TimeMidnight * 0.8f * landx);

#ifdef SKY_LIGHTING
color_sunlight         *= 9.8f * mix(1.0f, 0.15f, TimeMidnight * landx) * SUNLIGHT_POWER;
#else
color_sunlight         *= 90.8f * mix(1.0f, 0.15f, TimeMidnight * landx) * SUNLIGHT_POWER;
#endif
color_sunlight         *= mix(1.0f, 0.0f, rainx); //rain
float sunlight_desat    = dot(color_sunlight, vec3(1.0f));
color_sunlight          = mix(color_sunlight, rodcolor * sunlight_desat, TimeMidnight * 0.8f);

#ifdef SKY_LIGHTING
color_reflected  *= 1.00f * mix(1.0f, 0.00f, TimeMidnight * landx) + TimeNoon * 0.333;
#else
color_reflected  *= 0.550f * mix(1.0f, 0.00f, TimeMidnight * landx) + TimeNoon * 0.133;
#endif
float color_reflected_desat = dot(color_reflected, vec3(1.0f));
color_reflected  = mix(color_reflected, color_reflected_desat * rodcolor, TimeMidnight * 0.8f * landx);
color_reflected  *= mix(1.0f, 0.0f, rainx); //rain
color_reflected  = max(color_reflected, vec3(0.0f));

#ifdef SKY_LIGHTING
color_ambfill    *= 1.090f * mix(1.0f, 0.00f, TimeMidnight * landx) + (0.5f - TimeNoon);
#else
color_ambfill    *= 1.290f * mix(1.0f, 0.00f, TimeMidnight * landx) + (1.0f - TimeNoon);
#endif
float color_ambfill_desat = dot(color_ambfill, vec3(1.0f));
color_ambfill 		= mix(color_ambfill, rodcolor * color_ambfill_desat, TimeMidnight * 0.8f * landx);
color_ambfill    *= mix(1.0f, 0.8f, rainx);

color_torchlight *= 7.50f * TORCHLIGHT_POWER;
color_lightning  *= 0.50f;

Specular_lightmap *= 55.0f;

color_sky *= 3.0f + TimeNoon * 4.95f;

color_water_torch *= 7.50f * TORCHLIGHT_POWER * 0.0f;

color_water_sky *= 18.59f * 0.0f;
float color_water_sky_gray = dot(color_water_sky, vec3(1.0f));
color_water_sky = mix(color_water_sky, color_water_sky_gray * rodcolor, TimeMidnight * 0.8f);

color_water_sunlight *= 0.0f;

color_nolight *= 0.001f;
float nolight_desat = dot(color_nolight, vec3(1.0f));
color_nolight = mix(color_nolight, nolight_desat * rodcolor, 0.8f);


//Add all light elements together
color.rgb =   color_skylight 
			+ color_sunlight 
			+ color_reflected 
			+ color_torchlight 
			+ color_lightning 
			+ color_nolight 
			+ Specular_lightmap 
			+ color_sky 
			+ color_water_sky 
			+ color_water_torch 
			+ color_water_sunlight 
			+ color_ambfill;


//Godrays
float GRa = 0.0f;

#ifdef GODRAYS
	const float grna = 3300.0f;

	 GRa = addGodRays(0.0f, Texcoord2, noiseX2*grna, noiseY2*grna, noiseX3*grna, noiseY3*grna, noiseX4*grna)/2.0;
	 GRa = 1.0f - GRa;
	 //GRa += allwaves;
	 //GRa *= allwaves;
	 //GRa += iswater*0.25f;
#endif


//Render crepuscular rays
#ifdef CREPUSCULAR_RAYS
	float crepRays = 0.0f;
	if (rainx < 0.9f) {
	//complex rays with fog
	      crepRays = DrawCrepuscularRays(noiseX2, landx);
		  crepRays = pow(crepRays, 0.4545f);
		  crepRays *= 0.5f;
	} else {
		  crepRays = 0.09f;
	}
#else
float crepRays = 0.090f;
#endif


/*
//Calculate Fog
float fogDensity = 0.000025f * FOG_DENSITY;

float fogFactor = exp(getDepth(texcoord.st) * fogDensity) - 1.0f;
vec3  fogColor = pow(ambient_color, vec3(1.0f)) * 1.0f + (crepRays * pow(sunlight_color, vec3(2.3f)) * 2111.0f);

color.rgb = mix(color.rgb, fogColor * 50.0f, min(fogFactor, 1.0f));
*/


//Scale colors
color.rgb *= 0.015f;
color.b *= 1.0f;
color.rgb = mix(color.rgb * 2.9f, color.rgb * 1.0f, landx);

#ifdef PRESERVE_COLOR_RANGE
	color.rgb *= 0.5f;
#endif


#ifdef SUN_GLOW

	color.rgb += DoSunGlow(/*scale*/ 0.3f, /*color*/ sunlight_color, /*power*/ 0.3f, /*landx*/ landx) * mix(10.7f, 0.0f, TimeMidnight);
			

#endif

//color.rgb = vec3(crepRays);

//gamma fix
color.r = pow(color.r, gammafix);
color.g = pow(color.g, gammafix);
color.b = pow(color.b, gammafix);



//Pseudo HDR auto exposure
	/*
	//Eye sky factor
	float eyeSkylightFactor = eyeBrightnessSmooth.y / 16.0f;
	
		  eyeSkylightFactor = min(eyeSkylightFactor, 16.0f) / 16.0f;
		  eyeSkylightFactor = pow(eyeSkylightFactor, 3.0f);
		  eyeSkylightFactor *= mix(1.0f, NIGHT_EXPOSURE_BIAS, TimeMidnight);
		  
	//color.rgb /= eyeSkylightFactor * 1.375f + 0.4f;
	*/

/*
//color clip check
	if (color.r >= 1.0f) {
		color.r = 0.0f;
	}
	if (color.g >= 1.0f) {
		color.g = 0.0f;
	}
	if (color.b >= 1.0f) {
		color.b = 0.0f;
	}
*/	

	
//Attempt to hide LDR artifacts
	color.r = pow(color.r, 1.0f/BANDING_FIX_FACTOR);
	color.g = pow(color.g, 1.0f/BANDING_FIX_FACTOR);
	color.b = pow(color.b, 1.0f/BANDING_FIX_FACTOR);
	
	
	
//Handle materials and store water sunlight occlusion in material IDs

	float mats = texture2D(composite, texcoord.st).r;
		  mats *= 255.0f;
		  
		  shadingWater *= 8.0f;
		  
		  mats = max(mats, shadingWater * iswater);
		  
		  mats /= 255.0f;
	

	
	

	gl_FragData[0] = vec4(color.rgb, 1.0f);
	gl_FragData[1] = texture2D(gdepth, texcoord.st);
	gl_FragData[2] = vec4(normal * 0.5f + 0.5f, 1.0f);
	gl_FragData[3] = vec4(mats, specularity, GRa, 1.0f);
	gl_FragData[4] = vec4(texture2D(composite, texcoord.st).b, allwaves * 0.7f, texture2D(gdepth, texcoord.st).b, 1.0f);
	//gl_FragData[4] = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	//gl_FragData[5] = vec4(GRa, allwaves * 0.7, sky_lightmap, 1.0f);
	//gl_FragData[6] = texture2D(gaux3, texcoord.st);
	//gl_FragData[6] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

}
