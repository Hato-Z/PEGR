#version 120





//to increase shadow draw distance, edit SHADOWDISTANCE and SHADOWHPL below. Both should be equal. Needs decimal point.
//disabling is done by adding "//" to the beginning of a line.





//ADJUSTABLE VARIABLES

#define BLURFACTOR 3.5
#define SHADOW_DARKNESS 1.70   // 1.0 Is defualt darkness. 2.0 is black shadows. 0.0 is no shadows.
#define SHADOWDISTANCE 80.0 
#define SHADOW_FILTER

/* SHADOWRES:1024 */
/* SHADOWHPL:80.0 */

#define SSAO
#define SSAO_LUMINANCE 0.0				// At what luminance will SSAO's shadows become highlights.
#define SSAO_STRENGTH 0.725              // Too much strength causes white highlights on extruding edges and behind objects
#define SSAO_NOISE true					// Randomize SSAO sample gathering. With noise enabled and SSAO_LOOP set to 1, you will see higher performance at the cost of fuzzy dots in shaded areas.
#define SSAO_NOISE_AMP 0.025				// Multiplier of noise. Higher values mean SSAO takes random samples from a larger radius. Big performance hit at higher values.
#define SSAO_MAX_DEPTH 0.5				// View distance of SSAO
#define SSAO_SAMPLE_DELTA 0.15			// Radius of SSAO shadows. Higher values cause more performance hit.
#define CORRECTSHADOWCOLORS				// Colors sunlight and ambient light correctly according to real-life.
#define SHADOWOFFSET 0.4				// Shadow offset multiplier. Values that are too low will cause artifacts.
//#define FXAA							// FXAA shader. Broken, but you can give it a try if you want.
#define GODRAYS
#define GODRAYS_EXPOSURE 0.20
#define GODRAYS_SAMPLES 6
#define GODRAYS_DECAY 0.85
#define GODRAYS_DENSITY 0.80

//END OF ADJUSTABLE VARIABLES






uniform sampler2D gcolor;
uniform sampler2D gdepth;
uniform sampler2D shadow;
uniform sampler2D gaux1;
uniform sampler2D gaux2;

varying vec4 texcoord;
varying vec4 lmcoord;

uniform int worldTime;

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferModelViewInverse;
uniform mat4 shadowProjection;
uniform mat4 shadowModelView;
uniform vec3 sunPosition;

uniform float near;
uniform float far;
uniform float viewWidth;
uniform float viewHeight;
uniform float rainStrength;
uniform float aspectRatio;



// Standard depth function.
float getDepth(const in vec2 coord) {
    return 2.0f * near * far / (far + near - (2.0f * texture2D(gdepth, coord).x - 1.0f) * (far - near));
}

//Auxilliary variables
float land = texture2D(gaux1, texcoord.st).b;
float noblur = texture2D(gaux1, texcoord.st).r;
vec3 sunPos = sunPosition;
vec2 Texcoord2 = texcoord.st;

//Crossfading conditionals

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;
float landx = clamp(land-0.5f, 0.0f, 0.1f)/0.1f;

//Calculate Time of Day

	float timefract = worldTime;

	float TimeSunrise  = ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0) + (1.0 - (clamp(timefract, 0.0, 4000.0)/4000.0));
	float TimeNoon     = ((clamp(timefract, 0.0, 4000.0)) / 4000.0) - ((clamp(timefract, 8000.0, 12000.0) - 8000.0) / 4000.0);
	float TimeSunset   = ((clamp(timefract, 8000.0, 12000.0) - 8000.0) / 4000.0) - ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0);
	float TimeMidnight = ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0) - ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0);



#ifdef SSAO

// Alternate projected depth (used by SSAO, probably AA too)
float getProDepth(const in vec2 coord) {
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


vec2 rand(const in vec2 coord) { //generating noise/pattern texture for dithering
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
  float garea = 1.5f; //gauss bell width    
  float diff = (depth1 - depth2) * 50.0f; //depth difference (0-100)
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

  float coordw = texCoord.x + dw/(depth*1.0f + 0.0f);
  float coordh = texCoord.y + dh/(depth*1.0f + 0.0f);
  float coordw2 = texCoord.x - dw/(depth*1.0f + 0.0f);
  float coordh2 = texCoord.y - dh/(depth*1.0f + 0.0f);
  
  //attenuate samples in spherical manner
  float spherical = 1.0f - distance(texCoord.xy, vec2(dw, dh)) * 0.1f / (SSAO_NOISE_AMP + SSAO_SAMPLE_DELTA);

  if (coordw  < 1.0f && coordw  > 0.0f && coordh < 1.0f && coordh  > 0.0f){
    vec2 coord = vec2(coordw , coordh);
    vec2 coord2 = vec2(coordw2, coordh2);
    int zfar = 0;
    temp = compareDepths(depth, getProDepth(coord),zfar);
/*
    //DEPTH EXTRAPOLATION:
    if (zfar > 0){
      temp2 = compareDepths(getProDepth(coord2),depth,zfar);
      temp += (1.0f-temp)*temp2; 
    }
	*/
  }

  return temp;  
}  



float getSSAOFactor(float  noiseX1o, 
					float  noiseY1o, 
					float  noiseX2o, 
					float  noiseY2o, 
					float  noiseX3o, 
					float  noiseY3o, 
					float  noiseX4o, 
					float  noiseY4o) {

  float incx = 1.0f / viewWidth * SSAO_SAMPLE_DELTA;
  float incy = 1.0f / viewHeight * SSAO_SAMPLE_DELTA;
  
  noiseX1o *= SSAO_NOISE_AMP*50.0f;
  noiseY1o *= SSAO_NOISE_AMP*50.0f;
  noiseX2o *= SSAO_NOISE_AMP*50.0f;
  noiseY2o *= SSAO_NOISE_AMP*50.0f;
  noiseX3o *= SSAO_NOISE_AMP*50.0f;
  noiseY3o *= SSAO_NOISE_AMP*50.0f;
  noiseX4o *= SSAO_NOISE_AMP*50.0f;
  noiseY4o *= SSAO_NOISE_AMP*50.0f;
  
  
  
	vec2 noise1 = rand(texCoord)*4.0; 
	
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
  float aoMult = SSAO_STRENGTH/3.0f;
  float aaDiff = (1.0f + 2.0f / 1.0f); // 1.0 is samples

    float npw  = (pw + 0.09f * noiseX1o) / cdepth;
    float nph  = (ph + 0.09f * noiseY1o) / cdepth;
	

	float npw2  = (pw*2.0f + 0.09f * noiseX2o) / cdepth;
    float nph2  = (ph*2.0f + 0.09f * noiseY2o) / cdepth;
	
	float npw3  = (pw*3.0f + 0.09f * noiseX3o) / cdepth;
    float nph3  = (ph*3.0f + 0.09f * noiseY3o) / cdepth;
	
	float npw4  = (pw*4.0f + 0.09f * noiseX4o) / cdepth;
    float nph4  = (ph*4.0f + 0.09f * noiseY4o) / cdepth;

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
		
	
	ao += calAO(depth, npw4, nph4) * aoMult/2.5f;
    ao += calAO(depth, npw4, -nph4) * aoMult/2.5f;
    ao += calAO(depth, -npw4, nph4) * aoMult/2.5f;
    ao += calAO(depth, -npw4, -nph4) * aoMult/2.5f;
	/*

	ao += calAO(depth, 2.0*npw2, 2.0*nph2) * aoMult/3.0;
    ao += calAO(depth, 2.0*npw2, -2.0*nph2) * aoMult/3.0;
    ao += calAO(depth, -2.0*npw2, 2.0*nph2) * aoMult/3.0;
    ao += calAO(depth, -2.0*npw2, -2.0*nph2) * aoMult/3.0;
	
	ao += calAO(depth, 3.0*npw3, 3.0*nph3) * aoMult/3.5;
    ao += calAO(depth, 3.0*npw3, -3.0*nph3) * aoMult/3.5;
    ao += calAO(depth, -3.0*npw3, 3.0*nph3) * aoMult/3.5;
    ao += calAO(depth, -3.0*npw3, -3.0*nph3) * aoMult/3.5;
	
	 ao += calAO(depth, 4.0*npw4, 4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, 4.0*npw4, -4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, -4.0*npw4, 4.0*nph4) * aoMult/4.0;
    ao += calAO(depth, -4.0*npw4, -4.0*nph4) * aoMult/4.0;
	*/
	
	ao /= 3.0f;
	ao = 1.0f-ao;	
  ao = clamp(ao, 0.0f, 0.5f) * 2.0f;
	
  return ao;
}

#endif


#ifdef FXAA

vec4 fxaa() { 
	float pw = 1.0 / viewWidth;
	float ph = 1.0 / viewHeight;
	
	vec3 rgbNW = texture2D(gcolor, texcoord.xy + vec2(-pw,-ph)).xyz;
	vec3 rgbNE = texture2D(gcolor, texcoord.xy + vec2(pw,-ph)).xyz;
	vec3 rgbSW = texture2D(gcolor, texcoord.xy + vec2(-pw,ph)).xyz;
	vec3 rgbSE = texture2D(gcolor, texcoord.xy + vec2(pw,ph)).xyz;
	vec3 rgbM =  texture2D(gcolor, texcoord.xy).xyz;
	
	vec3 luma = vec3(0.299, 0.587, 0.114);
	float lumaNW = dot(rgbNW, luma);
	float lumaNE = dot(rgbNE, luma);
	float lumaSW = dot(rgbSW, luma);
	float lumaSE = dot(rgbSE, luma);
	float lumaM = dot(rgbM, luma);
	
	float lumaMin = min(lumaM, min(min(lumaNW,lumaNE),min(lumaSW,lumaSE)));
	float lumaMax = max(lumaM, max(max(lumaNW,lumaNE),max(lumaSW,lumaSE)));
	
	vec2 dir;
	dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
	dir.y = ((lumaNW + lumaSW) - (lumaNE + lumaSE));
	
	float dirReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 * 1.0/8.0), 1.0/128.0);
	
	float rcpDirMin = 1.0 /(min(abs(dir.x), abs(dir.y)) + dirReduce);
	
	dir = min(vec2(4.0,4.0),max(vec2(-4.0,-4.0),dir * rcpDirMin)) * vec2(pw,ph);
	
	vec3 rgbA = (1.0/2.0) * (texture2D(gcolor, texcoord.xy + dir * (1.0/3,0 - 0.5)).xyz + 
							texture2D(gcolor, texcoord.xy + dir * (2.0/3.0 - 0.5)).xyz);
							
	vec3 rgbB = rgbA * (1.0/2.0) + (1.0/4.0) * (texture2D(gcolor, texcoord.xy + dir * (0.0/3.0 - 0.5)).xyz + 
												texture2D(gcolor, texcoord.xy + dir * (3.0/3.0 - 0.5)).xyz);
	
	float lumaB = dot(rgbB, luma);
	
	if((lumaB < lumaMin) || (lumaB > lumaMax)) {
		return vec4(rgbA , 1.0);
	}
	return vec4(rgbB, 1.0);
}



#endif


#ifdef GODRAYS



	float addGodRays(in float nc, in vec2 tx, in float noise, in float noise2, in float noise3, in float noise4, in float noise5) {
			float GDTimeMult = 0.0f;
			if (sunPos.z > 0.0f) {
				sunPos.z = -sunPos.z;
				sunPos.x = -sunPos.x;
				sunPos.y = -sunPos.y;
				GDTimeMult = 0.0;
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
			
			for (int i = 0; i < 2; i++) {
			
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
			
			float bubble = distance(vec2(delta.x*aspectRatio, delta.y), vec2(0.0f, 0.0f))*1.0f;
				  bubble = clamp(bubble, 0.0f, 1.0f);
				  bubble = 1.0f - bubble;
				  
			return (nc + GODRAYS_EXPOSURE * (colorGD*bubble))*GDTimeMult;
        
	}
#endif 










void main() {

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



	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * texture2D(gdepth, texcoord.st).x - 1.0f, 1.0f);
	fragposition /= fragposition.w;
	
	#ifdef SHADOWDISTANCE
	float drawdistance = SHADOWDISTANCE;
	float drawdistancesquared = pow(drawdistance, 2);
	#endif
	
	float distance = sqrt(fragposition.x * fragposition.x + fragposition.y * fragposition.y + fragposition.z * fragposition.z);

	float shading = 1.0f;
	
	vec4 worldposition = vec4(0.0);

	if (distance < drawdistance && distance > 0.1f) {
		
		worldposition = gbufferModelViewInverse * fragposition;

		float xzDistanceSquared = worldposition.x * worldposition.x + worldposition.z * worldposition.z;
		float yDistanceSquared  = worldposition.y * worldposition.y;
		
		if (yDistanceSquared < drawdistancesquared) {
			
			worldposition = shadowModelView * worldposition;
			float comparedepth = -worldposition.z;
			worldposition = shadowProjection * worldposition;
			worldposition /= worldposition.w;
			
			worldposition.st = worldposition.st * 0.5f + 0.5f;
				
			if (comparedepth > 0.0f && worldposition.s < 1.0f && worldposition.s > 0.0f && worldposition.t < 1.0f && worldposition.t > 0.0f){
				float shadowMult = min(1.0f - xzDistanceSquared / drawdistancesquared, 1.0f) * min(1.0f - yDistanceSquared / drawdistancesquared, 1.0f);
				shadowMult = pow(shadowMult, 0.3f);
				float sampleDistance = 1.0f / 2048.0f;
					
					
					
			
					
				
				
					float zoffset = 0.00f;
					float offsetx = 0.0000f*BLURFACTOR*SHADOWOFFSET;
					float offsety = 0.0004f*BLURFACTOR*SHADOWOFFSET;
				
					
					//shadow filtering
					
					float step = 1.0f/2048.0f;
					float shadowdarkness = 0.5f*SHADOW_DARKNESS;
					float diffthresh = 0.75f;
					float bluramount = 0.00009f*BLURFACTOR;
					
					
					
					
					
					
					shading += -0.02 - shadowMult * (clamp(comparedepth - (0.05f + (texture2D(shadow, worldposition.st + vec2(1.0f*bluramount+offsetx, 1.0f*bluramount+offsety)).z) * (256.0f - 0.05f)) - zoffset, 0.0f, diffthresh)/diffthresh * shadowdarkness - zoffset);

				#ifdef SHADOW_FILTER
				
						shading = 1.0f;
			
					
						
							
							float xfadescale = 3.4f;
							
							/*
							if (rainStrength > 0.1) {
								xfadescale = 0.0;
							}
							*/
							
							xfadescale = mix(xfadescale, 5.0f, rainx);


							float trans = 0.0005f * xfadescale;
							
							float neg = 0.0f;
							
							float diffthresh2 = 0.7f;
							
							
					//bluramount *= xfadescale;
					
				
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*0.5, noiseY2*xfadescale + offsety + trans*0.5)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*9.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*0.5, noiseY3*xfadescale + offsety - trans*0.5)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*9.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*0.5, noiseY4*xfadescale + offsety + trans*0.5)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*9.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*0.5, noiseY5*xfadescale + offsety - trans*0.5)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*9.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans, noiseY2*xfadescale + offsety + trans)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*8.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans, noiseY3*xfadescale + offsety - trans)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*8.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans, noiseY4*xfadescale + offsety + trans)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*8.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans, noiseY5*xfadescale + offsety - trans)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*8.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*2.0, noiseY2*xfadescale + offsety + trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*2.0, noiseY3*xfadescale + offsety - trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*2.0, noiseY4*xfadescale + offsety + trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*2.0, noiseY5*xfadescale + offsety - trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*0.0, noiseY2*xfadescale + offsety + trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*0.0, noiseY3*xfadescale + offsety - trans*2.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx + trans*2.0, noiseY4*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*2.0, noiseY5*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*7.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*3.0, noiseY2*xfadescale + offsety + trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*3.0, noiseY3*xfadescale + offsety - trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*3.0, noiseY4*xfadescale + offsety + trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*3.0, noiseY5*xfadescale + offsety - trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*3.0, noiseY2*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx - trans*3.0, noiseY3*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*0.0, noiseY4*xfadescale + offsety + trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*0.0, noiseY5*xfadescale + offsety - trans*3.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*6.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*4.0, noiseY2*xfadescale + offsety + trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*4.0, noiseY3*xfadescale + offsety - trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*4.0, noiseY4*xfadescale + offsety + trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*4.0, noiseY5*xfadescale + offsety - trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*0.0, noiseY2*xfadescale + offsety + trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*0.0, noiseY3*xfadescale + offsety - trans*4.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx + trans*4.0, noiseY4*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*4.0, noiseY5*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*5.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*5.0, noiseY2*xfadescale + offsety + trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*5.0, noiseY3*xfadescale + offsety - trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*5.0, noiseY4*xfadescale + offsety + trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*5.0, noiseY5*xfadescale + offsety - trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
				
				
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*0.0, noiseY2*xfadescale + offsety + trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*0.0, noiseY3*xfadescale + offsety - trans*5.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx + trans*5.0, noiseY4*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*5.0, noiseY5*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*4.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*6.0, noiseY2*xfadescale + offsety + trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*6.0, noiseY3*xfadescale + offsety - trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*6.0, noiseY4*xfadescale + offsety + trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*6.0, noiseY5*xfadescale + offsety - trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*6.0, noiseY2*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx - trans*6.0, noiseY3*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*0.0, noiseY4*xfadescale + offsety + trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*0.0, noiseY5*xfadescale + offsety - trans*6.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*3.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*7.0, noiseY2*xfadescale + offsety + trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx + trans*7.0, noiseY3*xfadescale + offsety - trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*7.0, noiseY4*xfadescale + offsety + trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*7.0, noiseY5*xfadescale + offsety - trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX2*xfadescale + offsetx + trans*7.0, noiseY2*xfadescale + offsety + trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX3*xfadescale + offsetx - trans*7.0, noiseY3*xfadescale + offsety - trans*0.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX4*xfadescale + offsetx - trans*0.0, noiseY4*xfadescale + offsety + trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(noiseX5*xfadescale + offsetx - trans*0.0, noiseY5*xfadescale + offsety - trans*7.0)).z) * (256.0 - 0.05)) - zoffset, neg, diffthresh2)/diffthresh2 * shadowdarkness - zoffset)*2.0;
					

					shading = (shading/57.0)/5.0 + 0.95;
					
					shading *= 10.0; 
					shading += 1.7;
					
					float shadingsharp = 1.0;
					shadingsharp += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(offsetx, offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * (shadowdarkness*0.5) - zoffset);

				
					

					shadingsharp *= 0.38;
					shadingsharp -= 0.1;
					
					//shadingsharp = 0.4;
					
					
					/*
										if (rainStrength > 0.1) {
											
											shading = 2.2;
											shadingsharp = 0.8;
										}
										*/
										
										//raining
										shading = mix(shading, shading*0.4, rainx);
										shadingsharp = mix(shadingsharp, 0.7, rainx);
										
										shading *= shadingsharp;
										shading *= 0.44;
										shading = clamp(shading, 0.0, 1.0);
										
									
									
										
										
				#endif
					
					//remove shadows from sky
					shading = mix(1.0, shading, landx);
					
					
					
					
					
				
			
			}
		}
	}
	
//overbright shading
shading = mix(shading, shading*1.15 + 0.23, landx);
 
	vec4 color = texture2D(gcolor, texcoord.st);
	
	
	
//Determine what is being illuminated by the sun. 1 = sunlight. 0 = shadow.
float sun_amb = clamp(((shading-0.16)*3.0-1.1), 0.0, 1.0);
float sun_amb2 = clamp(((shading-0.21)*3.0-1.1), 0.0, 1.0);
	  sun_amb2 = mix(sun_amb2, 0.0f, rainx);

//raining
//sun_amb = mix(sun_amb, 0.4, rainx);
	
	

  color.rgb = mix(color.rgb, (color.rgb + 0.015), clamp(sun_amb * 2.0 - 1.0, 0.0, 1.0));
  
  
  
  
  float AO = 1.0;
#ifdef SSAO
	
  float lum = dot(color.rgb, vec3(1.0));
  vec3 luminance = vec3(lum);
  vec3 color_hold = vec3(color.rgb);
  AO *= getSSAOFactor(noiseX2, noiseY2, noiseX3, noiseY3, noiseX4, noiseY4, noiseX5, noiseY5);
  
  AO = mix(AO, AO * 0.5 + 0.5, sun_amb);
  
  if (land < 0.5) {
	AO = 1.0;
  }

  
  color.r *= AO;
  color.g *= AO;
  color.b *= AO;  

  
  color.r *= (AO*0.5 + 0.5);
  color.g *= (AO*0.5 + 0.5);
  color.b *= (AO*0.5 + 0.5);
  

  
  //color.rgb *= (clamp((AO*1.0), 0.0, 1.0));
  
#endif



 

#ifdef CORRECTSHADOWCOLORS

	vec3 sunrise_sun = vec3(1.0f, 0.73f, 0.40f) * TimeSunrise;
	
	vec3 sunrise_amb = vec3(0.9f, 0.9f, 0.9f) * TimeSunrise;
	
	vec3 noon_sun = vec3( mix(1.0f, 1.0f, rainx), mix(0.95f, 1.0f, rainx), mix(0.70f, 1.0f, rainx)) * TimeNoon;
	
	vec3 noon_amb = vec3( mix(0.80f, 1.0f, rainx), mix(0.93f, 1.0f, rainx), mix(1.00f, 1.0f, rainx)) * TimeNoon;
	
	vec3 sunset_sun = vec3(0.9f, 0.73f, 0.40f) * TimeSunset;
	
	vec3 sunset_amb = vec3(0.9f, 0.9f, 0.9f) * TimeSunset;
	
	vec3 midnight_sun = vec3(0.6f, 0.65f, 0.8f) * TimeMidnight;
	
	vec3 midnight_amb = vec3(0.8f, 0.85f, 1.0f) * TimeMidnight;

	vec3 sunlight = sunrise_sun + noon_sun + sunset_sun + midnight_sun;
	
	vec3 ambient = sunrise_amb + noon_amb + sunset_amb + midnight_amb;
	
	vec3 colorland = mix(color.rgb * ambient, color.rgb * sunlight, sun_amb);
	
	/*
	color.r = mix(color.r, colorland_r, landx);
	color.g = mix(color.g, colorland_g, landx);
	color.b = mix(color.b, colorland_b, landx);
	*/
	
	color.rgb = colorland;
	//color.r = mix(colorland_r, colorland_r, landx);
	//color.g = mix(colorland_g, colorland_g, landx);
	//color.b = mix(colorland_b, colorland_b, landx);

	//Sky saturation
	const float skysat_sunrise = 0.0f * TimeSunrise;
	const float skysat_noon = -0.3f * TimeNoon;
	const float skysat_sunset = 0.0f * TimeSunset;
	const float skysat_midnight = 0.0f * TimeMidnight;
	const float skysat = skysat_sunrise + skysat_noon + skysat_sunset + skysat_midnight;
	
	
		float colorskyclearr = (((color.r + (color.r * skysat) - ((color.g + color.b)*skysat)/2.0) * 1.5 - 0.3) * (TimeSunrise))   +   (( (color.r + (color.r * skysat) - ((color.g + color.b)*skysat)/2.0) * 2.25 - 1.0) * (TimeNoon))   +   (((color.r + (color.r * skysat) - ((color.g + color.b)*skysat)/2.0) * 1.5 - 0.3) * (TimeSunset))   +   (color.r * TimeMidnight);
		float colorskyclearg = (((color.g + (color.g * skysat) - ((color.r + color.b)*skysat)/2.0) * 1.5 - 0.3) * (TimeSunrise))   +   (( (color.g + (color.g * skysat) - ((color.r + color.b)*skysat)/2.0) * 2.15 - 1.0) * (TimeNoon))   +   (((color.g + (color.g * skysat) - ((color.r + color.b)*skysat)/2.0) * 1.5 - 0.3) * (TimeSunset))   +   (color.g * TimeMidnight);
		float colorskyclearb = (((color.b + (color.b * skysat) - ((color.g + color.r)*skysat)/2.0) * 1.7 - 0.3) * (TimeSunrise))   +   (( (color.b + (color.b * skysat) - ((color.g + color.r)*skysat)/2.0) * 2.25 - 1.0) * (TimeNoon))   +   (((color.b + (color.b * skysat) - ((color.g + color.r)*skysat)/2.0) * 1.7 - 0.3) * (TimeSunset))   +   (color.b * TimeMidnight);


			float colorskyrainr = ((color.r * 1.20 + 0.0) * (TimeSunrise))   +   ((color.r * 1.2 + 0.1) * (TimeNoon))   +   ((color.r * 1.2 + 0.0) * (TimeSunset))   +   (color.r * TimeMidnight);
			float colorskyraing = ((color.g * 1.20 + 0.0) * (TimeSunrise))   +   ((color.g * 1.2 + 0.1) * (TimeNoon))   +   ((color.g * 1.2 + 0.0) * (TimeSunset))   +   (color.g * TimeMidnight);
			float colorskyrainb = ((color.b * 1.20 + 0.0) * (TimeSunrise))   +   ((color.b * 1.2 + 0.1) * (TimeNoon))   +   ((color.b * 1.2 + 0.0) * (TimeSunset))   +   (color.b * TimeMidnight);
		
			float colorskyr = mix(colorskyclearr, colorskyrainr, rainx);
			float colorskyg = mix(colorskyclearg, colorskyraing, rainx);
			float colorskyb = mix(colorskyclearb, colorskyrainb, rainx);
			
			color.r = mix(colorskyr, color.r, landx);
			color.g = mix(colorskyg, color.g, landx);
			color.b = mix(colorskyb, color.b, landx);
			
	
	
/*
if (rainStrength > 0.1) {
	color.rgb *= 0.9;
}
*/




//raining
color.rgb = mix(color.rgb, color.rgb*0.8, rainx);
#endif

//brighter sky

color.rgb = mix(color.rgb*1.5, color.rgb, landx);


//PSEUDO SUN////////////////////

vec3 sP = sunPosition;

			vec2 lPos = sP.xy / -sP.z;
			lPos.y *= 1.2f;
			lPos.x *= 0.70f;
			lPos = (lPos + 1.0f)/2.0f;
			
			
			
			//Rain glow
			const vec2 glow1scale = vec2(0.5f, 0.5f);
			const float glow1pow = 1.5f;
			const float glow1fill = 1.0f;
			const float glow1offset = -2.0f;
			vec2 glow1pos = vec2(  ((1.0 - lPos.x)*(glow1offset + 1.0) - (glow1offset*0.5))  *aspectRatio*glow1scale.x,  ((1.0 - lPos.y)*(glow1offset + 1.0) - (glow1offset*0.5))  *glow1scale.y);
			
			
			float glow1 = distance(glow1pos, vec2(texcoord.s*aspectRatio*glow1scale.x, texcoord.t*glow1scale.y));
				  glow1 = 0.5 - glow1;
				  glow1 = clamp(glow1*glow1fill, 0.0, 1.0);
				  glow1 = sin(glow1*1.57075);
				  glow1 *= (1.0f - landx);
				  glow1 = pow(glow1, 2.5f);
				  
				  glow1 *= glow1pow;
				  glow1 *= rainx;
				  
				  	color.r += glow1*sunlight.r * (1.0f - (TimeMidnight*0.7f));
					color.g += glow1*sunlight.g * (1.0f - (TimeMidnight*0.7f));
					color.b += glow1*sunlight.b * (1.0f - (TimeMidnight*0.7f));				
			
			
			
			//Sunny glow
			const vec2 glow2scale = vec2(0.8f, 0.8f);
			const float glow2pow = 1.29f;
			const float glow2fill = 1.0f;
			const float glow2offset = -2.0f;
			vec2 glow2pos = vec2(  ((1.0 - lPos.x)*(glow2offset + 1.0) - (glow2offset*0.5))  *aspectRatio*glow2scale.x,  ((1.0 - lPos.y)*(glow2offset + 1.0) - (glow2offset*0.5))  *glow2scale.y);
			
			
			float glow2 = distance(glow2pos, vec2(texcoord.s*aspectRatio*glow2scale.x, texcoord.t*glow2scale.y));
				  glow2 = 0.5 - glow2;
				  glow2 = clamp(glow2*glow2fill, 0.0, 1.0);
				  glow2 = sin(glow2*1.57075);
				  glow2 *= (1.0f - landx);
				  glow2 = pow(glow2, 2.5f);
				  
				  glow2 *= glow2pow;
				  glow2 *= 1.0 - rainx;
				  
				  	color.r += glow2*sunlight.r * (1.0f - (TimeMidnight*0.7f));
					color.g += glow2*sunlight.g * (1.0f - (TimeMidnight*0.7f));
					color.b += glow2*sunlight.b * (1.0f - (TimeMidnight*0.7f));	

///////////////////////////////


//Apply Shadows
color.rgb *= shading;
color.rgb *= 0.9f;

float GRa = 0.0f;




#ifdef GODRAYS
	const float grna = 3300.0f;

		  GRa = addGodRays(0.0f, Texcoord2, noiseX3*grna, noiseX4*grna, noiseY4*grna, noiseX2*grna, noiseY2*grna)/2.0;
		  GRa *= 1.0 - rainx;
#endif



//Adjust blue channel based on weather
color.b *= 1.10;

color.b = mix(color.b, color.b*1.2, rainx);


//boost red channel
color.r *= 1.05;




//Average luminosity for HDR
float HDRlum = (color.r + color.g + color.b)/3.0f;
	  HDRlum *= shading;
	  HDRlum = shading;
	  



	//gl_FragData[0] = texture2D(gcolor, texcoord.st);
	gl_FragData[1] = texture2D(gdepth, texcoord.st);
	gl_FragData[3] = vec4(color.rgb, 1.0);
	gl_FragData[4] = vec4(noblur, 1.0 - GRa, land, 0.0f);
}

