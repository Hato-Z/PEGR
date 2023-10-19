#version 120



//to increase shadow draw distance, edit SHADOWDISTANCE and SHADOWHPL below. Both should be equal. Needs decimal point.
//disabling is done by adding "//" to the beginning of a line.



//ADJUSTABLE VARIABLES

#define BLURFACTOR 3.5
#define SHADOW_DARKNESS 0.6
#define SUNLIGHTAMOUNT 1.1
#define SHADOWDISTANCE 75.0 
#define SHADOW_CLAMP 0.5
#define SHADOW_RES 1024.0

/* SHADOWRES:1024 */
/* SHADOWHPL:75.0 */

  #define SSAO
  #define SSAO_LUMINANCE 0.0				// At what luminance will SSAO's shadows become highlights.
  #define SSAO_STRENGTH 1.75               // Too much strength causes white highlights on extruding edges and behind objects
  #define SSAO_LOOP 1						// Integer affecting samples that are taken to calculate SSAO. Higher values mean more accurate shadowing but bigger performance impact
  #define SSAO_MAX_DEPTH 0.75				// View distance of SSAO
  #define SSAO_SAMPLE_DELTA 0.4			// Radius of SSAO shadows. Higher values cause more performance hit.
  #define CORRECTSHADOWCOLORS				// Colors sunlight and ambient light correctly according to real-life. 
  #define SHADOWOFFSET 0.0				// Shadow offset multiplier. Values that are too low will cause artefacts.
  #define GODRAYS
  #define GODRAYS_EXPOSURE 0.10
  #define GODRAYS_SAMPLES 6
  #define GODRAYS_DECAY 0.95
  #define GODRAYS_DENSITY 0.65

uniform int fogMode;
  
  #define BUMPMAPPWR 0.55			//fine for x64 texturepack, better if you lower this when using lower res texturepacks
  //#define CLAY_RENDER

//#define PRESERVE_COLOR_RANGE


//END OF ADJUSTABLE VARIABLES






uniform sampler2D gcolor;
uniform sampler2D depthtex0;
uniform sampler2D gnormal;
uniform sampler2D shadow;
uniform sampler2D gaux1;

varying vec4 texcoord;
varying vec4 lmcoord;
varying vec3 lightVector;

uniform int worldTime;

uniform mat4 gbufferProjection;
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


//Calculate Time of Day

varying float TimeMidnight;
varying float TimeSunset;
varying float TimeNoon;
varying float TimeSunrise;

//colors
varying vec3 sunlight_color;
varying vec3 ambient_color;
varying vec3 skycolor;

varying vec3 sunlight;

float edepth(vec2 coord) {
return texture2D(depthtex0,coord).z;
}
float luma(vec3 color) {
return dot(color.rgb,vec3(0.299, 0.587, 0.114));
}
float ld(float depth) {
    return (2.0 * near) / (far + near - depth * (far - near));
}

vec2 texel = vec2(1.0/viewWidth,1.0/viewHeight);
//Auxilliary variables
vec3 aux = texture2D(gaux1, texcoord.st).rgb;

float	land 			 = aux.b;
//float	noblur 			 = texture2D(gaux1, texcoord.st).r;
vec3	sunPos			 = sunPosition;
vec2 	Texcoord2		 = texcoord.st;
float 	iswater			 = 0.0;
vec3 	normal         	 = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
float 	translucent		 = 0.0;
float  pixeldepth = texture2D(depthtex0,texcoord.xy).x;
float totalspec = aux.r*3.0;

//Crossfading conditionals

float rainx = clamp(rainStrength, 0.0f, 1.0f);
float wetx  = clamp(wetness, 0.0f, 1.0f);

float pw = 1.0/ viewWidth;
float ph = 1.0/ viewHeight;
//Lightmaps

float sky_lightmap = 1.0;
float torch_lightmap = pow(aux.b,1.5);
float lightning_lightmap = 0.0;


float shadowexit = 0.0;


			


#ifdef SSAO

// Alternate projected depth (used by SSAO, probably AA too)
float getProDepth(vec2 coord) {
	float depth = texture2D(depthtex0, coord).x;
	return ( 2.0f * near ) / ( far + near - depth * ( far - near ) );
}

float znear = near; //Z-near
float zfar = far; //Z-far

float diffarea = 0.6f; //self-shadowing reduction
float gdisplace = 0.35f; //gauss bell center

//bool noise = SSAO_NOISE; //use noise instead of pattern for sample dithering?
bool onlyAO = false; //use only ambient occlusion pass?

vec2 texCoord = texcoord.st;



float compareDepths(in float depth1, in float depth2) {  
  float garea = 8.5f; //gauss bell width    
  float diff = (depth1 - depth2) * 100.0f; //depth difference (0-100)
  //reduce left bell width to avoid self-shadowing 
  
  if (diff < gdisplace) {
    garea = diffarea;
  } 


  float gauss = pow(2.7182f,-2.0f*(diff-gdisplace)*(diff-gdisplace)/(garea*garea));
  return gauss;
} 

float calAO(float depth, vec2 coord) {  
  float temp = 0.0f;
  vec2 coord2 = texcoord.xy + coord/clamp(0.15+depth*0.2,0.0,0.22);

    temp = compareDepths(depth, getProDepth(coord2));


  return temp;  
}  



float getSSAOFactor() {

  vec2 inc = texel;
  
	float depth = ld(pixeldepth);
	
  if (depth > SSAO_MAX_DEPTH) {
    return 1.0f;
  }
  float cdepth = pixeldepth;

	float ao = 0.0;
	float s;
	


  float aoMult = SSAO_STRENGTH;


    vec2 np = inc;
	
	float noiseX4 = clamp(fract(sin(dot(Texcoord2 ,vec2(16.9898f,38.633f))) * 41178.5453f),0.0f,1.0f)*2.0f-1.0f;
						

	

    ao += calAO(depth, np) * aoMult;
    ao += calAO(depth, vec2(np.x,-np.y)) * aoMult;
    ao += calAO(depth, vec2(-np.x,np.y)) * aoMult;
    ao += calAO(depth, -np) * aoMult;
	
	ao += calAO(depth, np*2.0) * aoMult/1.5;
    ao += calAO(depth, vec2(np.x,-np.y)*2.0) * aoMult/1.5;
    ao += calAO(depth, vec2(-np.x,np.y)*2.0) * aoMult/1.5;
    ao += calAO(depth, -np*2.0) * aoMult/1.5;
	
    ao += calAO(depth, np*3.0) * aoMult/2.0;
    ao += calAO(depth, vec2(np.x,-np.y)*3.0) * aoMult/2.0;
    ao += calAO(depth, vec2(-np.x,np.y)*3.0) * aoMult/2.0;
    ao += calAO(depth, -np*3.0) * aoMult/2.0;
	
	
	
	ao /= 16.0f;
	ao = 1.0f-ao;	
  ao = clamp(ao, 0.0f, 0.5f) * 2.0f;
	
  return ao;
}

#endif


#ifdef GODRAYS



	float addGodRays(vec2 lightPos, in float nc, in vec2 tx, in float noise, in float noise2, in float noise3, in float noise4, in float noise5) {
			float GDTimeMult = 0.0f;
			if (sunPos.z > 0.0f) {
				sunPos.z = -sunPos.z;
				sunPos.x = -sunPos.x;
				sunPos.y = -sunPos.y;
				GDTimeMult = TimeMidnight;
			} else {
				GDTimeMult = TimeSunrise + TimeNoon + TimeSunset;
			}
			//vec2 coord = tx;
			vec2 delta = (tx - lightPos) * GODRAYS_DENSITY / 2.0;
			delta *= -sunPos.z*0.01f;
			//delta *= -sunPos.z*0.01;
			float decay = -sunPos.z / 100.0f;
				 // decay *= -sunPos.z*0.01;
			float colorGD = 0.0f;
			for (int i = 0; i<2;i++) {
			

				
			
			
				tx -= delta;
				float sample = 0.0f;

					sample = step(texture2D(gaux1, clamp(tx + delta*noise,0.000001,0.999999)).g,0.01);
					sample += step(texture2D(gaux1, clamp(tx + delta*noise2,0.000001,0.999999)).g,0.01);
					sample += step(texture2D(gaux1, clamp(tx + delta*noise3,0.000001,0.999999)).g,0.01);
					sample += step(texture2D(gaux1, clamp(tx + delta*noise4,0.000001,0.999999)).g,0.01);
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


//
float callwaves(vec3 pos) {
float wsize = 9.0;
float wspeed = 0.25f;

float rs0 = abs(sin((worldTime*wspeed/5.0) + (pos.s*wsize) * 20.0)+0.2);
float rs1 = abs(sin((worldTime*wspeed/7.0) + (pos.t*wsize) * 27.0));
float rs2 = abs(sin((worldTime*wspeed/2.0) + (pos.t*wsize) * 60.0 - sin(pos.s*wsize) * 13.0)+0.4);
float rs3 = abs(sin((worldTime*wspeed/1.0) - (pos.s*wsize) * 20.0 + cos(pos.t*wsize) * 83.0)+0.1);

float wsize2 = 7.0;
float wspeed2 = 0.17f;

float rs0a = abs(sin((worldTime*wspeed2/4.0) + (pos.s*wsize2) * 24.0));
float rs1a = abs(sin((worldTime*wspeed2/11.0) + (pos.t*wsize2) * 77.0 )+0.3);
float rs2a = abs(sin((worldTime*wspeed2/6.0) + (pos.s*wsize2) * 50.0 - (pos.t*wsize2) * 23.0)+0.12);
float rs3a = abs(sin((worldTime*wspeed2/14.0) - (pos.t*wsize2) * 4.0 + (pos.s*wsize2) * 98.0));

float wsize3 = 3.0;
float wspeed3 = 0.3f;

float rs0b = abs(sin((worldTime*wspeed3/4.0) + (pos.s*wsize3) * 14.0));
float rs1b = abs(sin((worldTime*wspeed3/11.0) + (pos.t*wsize3) * 37.0));
float rs2b = abs(sin((worldTime*wspeed3/6.0) + (pos.t*wsize3) * 47.0 - cos(pos.s*wsize3) * 33.0 + rs0a + rs0b));
float rs3b = abs(sin((worldTime*wspeed3/14.0) - (pos.s*wsize3) * 13.0 + sin(pos.t*wsize3) * 98.0 + rs0 + rs1));

float waves = (rs1 * rs0 + rs2 * rs3)/2.0f;
float waves2 = (rs0a * rs1a + rs2a * rs3a)/2.0f;
float waves3 = (rs0b + rs1b + rs2b + rs3b)*0.25;


return (waves + waves2 + waves3)/3.0f;
}

vec2 rand(vec2 coord) { //generating noise/pattern texture for dithering
  const float width = 1.0f;
  const float height = 1.0f;
  float noiseX = ((fract(1.0f-coord.s*(width/2.0f))*0.25f)+(fract(coord.t*(height/2.0f))*0.75f))*2.0f-1.0f;
  float noiseY = ((fract(1.0f-coord.s*(width/2.0f))*0.75f)+(fract(coord.t*(height/2.0f))*0.25f))*2.0f-1.0f;

  //generate SSAO noise
  noiseX = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f))) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  noiseY = clamp(fract(sin(dot(coord ,vec2(12.9898f,78.233f)*2.0f)) * 43758.5453f),0.0f,1.0f)*2.0f-1.0f;
  
  return vec2(noiseX,noiseY);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void main() {
vec4 tpos = vec4(sunPos,1.0)*gbufferProjection;
			tpos = vec4(tpos.xyz/tpos.w,1.0);
			vec2 lightPos = tpos.xy/tpos.z;
			lightPos = (lightPos + 1.0f)/2.0f;
			
float rainmask = 0.0;
if(aux.g > 0.1 && aux.g < 0.3){
shadowexit = 1.0;
land = 1.0;
iswater = 0.0;
}

if(aux.g < 0.01) {
land = 0.0;
shadowexit = 1.0;
iswater =0.0;
}

if(aux.g > 0.01 && aux.g < 0.07) {
iswater = 1.0;
land = 1.0;
shadowexit = 0.0;
}

if(aux.g > 0.3 && aux.g < 0.5) {
iswater = 0.0;
land = 1.0;
shadowexit = 0.0;
translucent = 1.0;
}

if(aux.g > 0.8) {
iswater = 0.0;
land = 1.0;
shadowexit = 0.0;
}

if(aux.g > 0.6 && aux.g < 0.8) {
rainmask = 1.0;
iswater = 0.0;
land = 0.0;
shadowexit = 0.0;
}

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



	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * pixeldepth - 1.0f, 1.0f);
	fragposition /= fragposition.w;
	
	#ifdef SHADOWDISTANCE
	float drawdistance = SHADOWDISTANCE;
	float drawdistancesquared = pow(drawdistance, 2.0f);
	#endif
	
	float dist = length(fragposition.xyz);

	float shading = 1.0f;
	float shadingsharp = 1.0f;
	float diffthresh = SHADOW_CLAMP;
	
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
			
			float sample = 0.0;
			float isoccluded = 0.0;


	  int vpsize = 0;
	 
float isshadow = 0.0;
float ssample;
		float tfade = abs(clamp(worldTime*1.0,3500.0,8500.0)-5500.0)/2000.0;
		if (dist < drawdistance*1.75) {
		
			

				
			if (comparedepth > 0.f && worldposition.s < 0.98 && worldposition.s > 0.02 && worldposition.t < 0.98 && worldposition.t > 0.02){
			float shadowMult = 1.0;
					if (shadowexit > 0.1) {
					
					shading = 0.0;
					}
					else
					{
					
					//shadow filtering
					float vpsdepth = 0.0f;
						float vpsconfusion = 1.5;
					float zoffset = 0.00f;
					float offsetx = -0.0000f*BLURFACTOR*SHADOWOFFSET*(TimeSunset * 2.0f - 1.0f);
					float offsety = 0.0000f*BLURFACTOR*SHADOWOFFSET;
						
							
							float maxCompareDepth = 60.0f;
							vec2 vpsSpread = vec2(1.0/SHADOW_RES) * vpsconfusion;
							vec2 vpsc = vec2(-2.0f);
							float vpssamp = 0.0f;
							float vpscurve = 0.45f;
							
							for(int i = 0; i < 3; i++){
								vpsc.x = -2.0f;
									
								for(int j = 0; j < 3; j++){
									ssample = comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vpsSpread * vpsc).z) * (256.0 - 0.05));
									vpsdepth += pow((clamp(ssample, 0.0, maxCompareDepth)/(maxCompareDepth)), 1.0f/vpscurve);
									isshadow += clamp(ssample, 0.0, diffthresh)/(diffthresh);
									vpsc.x += 2.0f;
									vpssamp += 1.0f;
								}
								vpsc.y += 2.0f;
							}
							
							vpsdepth = pow(vpsdepth, vpscurve);
							vpsdepth /= vpssamp;
							vpsdepth *= 0.5;
							
						
						
						diffthresh = diffthresh * vpsdepth * 45.0f + diffthresh;
					
						float sfx = 7.5 * vpsdepth/SHADOW_RES;
						
						float sfxc = -4.0f;
						float sfyc = -4.0f;
						
						float sfsamp = 1.0f;
						
					
						
						//don't filter if not occluded/fully occluded
								if (isshadow > 0.9 && isshadow < 7.9 && vpsdepth > 0.007) {
								for(int i = 0; i < 9; i++){
								for(int j = 0; j < 9; j++){
									shadingsharp += (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(sfxc, sfyc)*sfx).z) * (256.0 - 0.05)), 0.0, diffthresh)/(diffthresh) - zoffset);
									sfxc += 1.0f;
									sfsamp += 1.0f;
								}
								sfxc = -4.0;
								sfyc += 1.0;
							}
							shadingsharp /= sfsamp;
							shading = 1.0-shadingsharp;
							isshadow = 1.0;
						}
						else {
						shading = 1.0-clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st).z) * (256.0 - 0.05)), 0.0,diffthresh)/diffthresh;
						isshadow = 0.0;
						}
						
						/*
						//add random samples
						float random = rand(texcoord.xy);
						//float random2 = randomize(worldpositionraw.xy);
						for (int i = 0; i < 5; i++) {						
						shadingsharp += (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st  + vec2(random)*sfx*(i-2)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/(diffthresh) - zoffset);
						sfsamp += 1.0f;
						}
						*/
						
						
						
}
	shading = clamp(shading + (1.0-tfade)*0.4,0.0,1.0);						
}
}

			

	
	
		  


	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////
	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////
	//////////////////DIRECTIONAL LIGHTING WITH NORMAL MAPPING////////////////////
float AO=1.0;

#ifdef SSAO
	if (land > 0.9 && iswater < 0.9) {

  AO = getSSAOFactor();
 
  
  AO = max(AO, 0.0f);
  }

#endif

			
			

float sunlight_direct = 1.0;
float sunlight_reflected = 0.0;
float spec = 0.0;
float sky_spec = 0.0;
float fresnel = 0.0;
vec3 specular;
vec3 npos;
if (land > 0.9 && iswater < 0.9) {
					npos = normalize(fragposition.xyz);

					
					specular = reflect(npos, normal);
					//simulating 3 (lights spots for sky fresnel fake
					 fresnel = (max(dot(specular, normalize(vec3(1.0,0.0,-1.0))),0.0)*0.5+max(dot(specular, normalize(vec3(-1.0,0.0,-1.0))),0.0))*0.6;			
						fresnel = clamp(pow(fresnel, 3.0f),0.0,0.8);
					
					float sunlight = max(dot(normal, lightVector),0.0);
	

					float direct  = sunlight;
						  				
						   if (dist < 60.0 && totalspec > 0.01) {	
					  
						  spec = max(dot(specular, lightVector), 0.0f);
						  spec = pow(spec, 20.0f);
						  spec = mix(0.0f, spec, clamp(shading, 0.0, 1.0));
						  spec *= 3.0;
						  spec *= totalspec;
						  sky_spec = fresnel*totalspec*0.5;
						  }

						  sunlight_direct = mix(1.0,direct,BUMPMAPPWR);
						  sunlight_direct = mix(sunlight_direct, 0.9, translucent);		//1.0 is looking weird on transparent stuff
						  
						  
					
					shading *= sunlight_direct;
				}
				
				
				shading = mix(shading,0.0,rainx);
 
//Albedo
vec3 color = texture2D(gcolor, texcoord.st).rgb;


vec3 torchcolor = vec3(1.0,0.675,0.415);


vec3 Specular_lightmap = clamp(spec * sunlight_color * shading * (1.0-rainx) + sky_spec*pow(sky_lightmap,4.0)*skycolor,0.0,1.1);
	 Specular_lightmap *= sky_lightmap;

vec3 Sunlight_lightmap = sunlight_color*shading*sky_lightmap*SUNLIGHTAMOUNT + ambient_color*(1.0-shading)*(1.0-SHADOW_DARKNESS)*sky_lightmap;
	 
	 
	 
vec3 Torchlight_lightmap = torch_lightmap *  torchcolor;

	 
vec3 LightningFlash_lightmap = vec3(lightning_lightmap *  0.8f, lightning_lightmap *  0.7f, lightning_lightmap *  1.0f);


//RAINWET
			float dampmask = clamp(sky_lightmap * 4.0f - 1.0f, 0.0f, 1.0f) * land * wetx;
			
			
			color = pow(color, vec3(mix(1.0f, 1.35f, dampmask)));

			
			
			
		

const float rspread = 0.30f;

float wave = 0.0;
if (iswater > 0.9) {
wave = callwaves(worldposition.xyz)*2.0-1.0;

float angle = dot(normalize(fragposition.xyz),normal);
wave = wave*(0.5+abs(angle)*0.5);


const float wnormalclamp = 0.05f;

float rdepth = pixeldepth;
		float waves = wave;
		float wnormal_x1 = texture2D(depthtex0, texcoord.st + vec2(pw, 0.0f)).x - texture2D(depthtex0, texcoord.st).x;
		float wnormal_x2 = texture2D(depthtex0, texcoord.st).x - texture2D(depthtex0, texcoord.st + vec2(-pw, 0.0f)).x;			
		float wnormal_x = 0.0f;
		
		if(abs(wnormal_x1) > abs(wnormal_x2)){
			wnormal_x = wnormal_x2;
		} else {
			wnormal_x = wnormal_x1;
		}
		wnormal_x /= 1.0f - rdepth;	

		wnormal_x = clamp(wnormal_x, -wnormalclamp, wnormalclamp);
		
		wnormal_x *= rspread;
		

			  
			  
		float wnormal_y1 = texture2D(depthtex0, texcoord.st + vec2(0.0f, ph)).x - texture2D(depthtex0, texcoord.st).x;
		float wnormal_y2 = texture2D(depthtex0, texcoord.st).x - texture2D(depthtex0, texcoord.st + vec2(0.0f, -ph)).x;		
		float wnormal_y;
		
		if(abs(wnormal_y1) > abs(wnormal_y2)){
			wnormal_y = wnormal_y2;
		} else {
			wnormal_y = wnormal_y1;
		}	
		wnormal_y /= 1.0f - rdepth;			

		wnormal_y = clamp(wnormal_y, -wnormalclamp, wnormalclamp);
		
		wnormal_y *= rspread;
		

//Calculate distance of objects behind water
float refractdist = 0.2 * 10.0f;

//Perform refraction
float refractamount = 500.1154f*0.35f*refractdist;
float refractamount2 = 0.0214f*0.05f*refractdist;
float refractamount3 = 0.214f*0.15f*refractdist;
float waberration = 0.105;

	vec3 refracted = vec3(0.0f);
	float refractedmask = 0.0;
	float bigWaveRefract = 600.0f * (1.0f - rdepth)*refractdist;
	float bigWaveRefractScale = 1000.0f * (1.0f - rdepth)*refractdist;
	
	vec2 bigRefract = vec2(wnormal_x*bigWaveRefract, wnormal_y*bigWaveRefract);
	
	vec2 refractcoord_r = texcoord.st;
	vec2 refractcoord_g = texcoord.st;
	vec2 refractcoord_b = texcoord.st;
	
	for (int i = 0; i < 1; ++i) {
			
	
			 refractcoord_r = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2( waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x,  waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y) * (waberration * 2.0f + 1.0f);

				
				refractcoord_r = refractcoord_r * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);

				
			
			refractcoord_r.s = clamp(refractcoord_r.s, 0.001f, 0.999f);
			refractcoord_r.t = clamp(refractcoord_r.t, 0.001f, 0.999f);	
			
			
			
			if (refractcoord_r.s > 1.0 || refractcoord_r.s < 0.0 || refractcoord_r.t > 1.0 || refractcoord_r.t < 0.0) {
					break;
				}
				
			

			
			refracted.rgb = texture2D(gcolor, refractcoord_r).rgb;
			
			
			refractedmask = texture2D(gaux1, refractcoord_r).g;
if(refractedmask > 0.01 && refractedmask < 0.07) {
refractedmask = 1.0;
}
else refractedmask = 0.0;
	
			}
			
	color.rgb = mix(color.rgb, refracted.rgb, vec3(refractedmask));
							
							
						  npos = normalize(fragposition.xyz);
						  specular = reflect(npos, normal+waves*0.025);
						  spec = max(dot(specular, lightVector), 0.0f);
						  spec = clamp(pow(spec, 30.0f),0.0,1.0)*(1.0+TimeMidnight);
						  spec = mix(0.0f, spec, clamp(shading, 0.0, 1.0));
						  
						  color.rgb += spec*sunlight_color;
						 
	
}




//Apply different lightmaps to image
if (land > 0.9 ) {
vec3 color_sunlight = color * Sunlight_lightmap;
vec3 color_torchlight = color * Torchlight_lightmap;

color_sunlight *= mix(1.0,0.6,TimeMidnight);


//Add all light elements together
color = color_sunlight + color_torchlight + Specular_lightmap ;

}
else 
{
color.rgb = mix(color.rgb,skycolor,rainx)*0.75;
}

//Godrays
float GRa = 0.0f;

#ifdef GODRAYS
	const float grna = 3300.0f;
			if (lightPos.x < 1.0 && lightPos.x > 0.0 && lightPos.y < 1.0 && lightPos.y > 0.0) {
			
	 GRa = addGodRays(lightPos,0.0f, Texcoord2, noiseX3*grna, noiseX4*grna, noiseY4*grna, noiseX2*grna, noiseY2*grna)/2.0;
}
#endif
color.rgb *= AO;
if (fogMode != 0) {
float fogfactor = pow(max(length(fragposition.xyz)/far-0.2+rainx*0.15,0.0),2.0-rainx*1.5);
color.rgb = mix(color.rgb,(sunlight+vec3(0.9-TimeMidnight*0.8))*0.5*(1.0-rainx*0.5),fogfactor*land);
}
vec3 rainclr = vec3(0.15,0.15,0.4);
color.rgb = mix(color.rgb,rainclr,rainmask*0.1);

float visiblesun = 0.0;
float temp;
int nb = 0;

//calculate sun occlusion (only on one pixel) 
if (texcoord.x < pw && texcoord.x < ph && lightPos.x < 1.0 && lightPos.x > 0.0 && lightPos.y < 1.0 && lightPos.y > 0.0) {
	for (int i = 0; i < 16;i++) {
		for (int j = 0; j < 16 ;j++) {
		temp = texture2D(gaux1,lightPos + vec2(pw*(i-8)*7.0,ph*(j-8)*7.0)).g;
		if (temp > 0.04) visiblesun += 0.0;
		else visiblesun += 1.0;
		nb += 1;
		}
	}
	visiblesun /= nb;

}

/* DRAWBUFFERS:NNN3N5 */
//color.rgb = vec3(isshadow);
    gl_FragData[5] = vec4(iswater,wave*0.5+0.5,GRa,visiblesun);
	gl_FragData[3] = vec4(color, land);
}
