#version 120


// Place two forward Slashes in front of the following '#define' lines in order to disable an option.
// MOTIONBLUR, HDR, and BOKEH_DOF are very beta shaders. Use at risk of weird results.
// MOTIONBLUR and BOKEH_DOF are not compatable with eachother. Shaders break when you enable both.
// GLARE is still a work in progress.
// BLOOM is currently broken.




  //#define BOKEH_DOF							//Cannot be applied to water
	//#define HQ_DOF								//Enable for higher quality DOF
	//#define TILT_SHIFT							//Tilt shift effect. Not meant for gameplay. Google "tilt shift" for more info.
	//#define TILT_SHIFT_SCALE 0.5				//Size of aperture. Higher values gives illusion of smaller world
  
  //#define CREPUSCULAR_RAYS
    //#define RAY_BRIGHTNESS 1.5
    
  #define FOG_DENSITY 0.35
  
  #define GODRAYS
	#define GODRAYS_EXPOSURE 0.02
	#define GODRAYS_SAMPLES 6
	#define GODRAYS_DECAY 0.99
	#define GODRAYS_DENSITY 0.30
	
  #define LENS								
	#define LENS_POWER 0.26
	
  //#define GLARE
	#define GLARE_AMOUNT 0.00006
	#define GLARE_RANGE 1.0
	#define GLARE2							//second pass of glare shader. More realistic light scattering.
  
  //#define CEL_SHADING
    #define CEL_SHADING_THRESHOLD 0.5
    #define CEL_SHADING_THICKNESS 0.002

  #define VIGNETTE
	#define VIGNETTE_STRENGTH 1.30
  
  //#define LOWLIGHT_EYE
  
  #define TONEMAP2	
  
  #define NIGHT_EXPOSURE_BIAS 0.1f
  //#define TONEMAP	
	//#define TONEMAP_FILMIC						
	//#define TONEMAP_COLOR_FILTER
  
  #define BRIGHTMULT 1.00               	// 1.0 = default brightness. Higher values mean brighter. 0 would be black.
  #define DARKMULT 0.08					// 0.0 = normal image. Higher values will darken dark colors.
  #define COLOR_BOOST	0.06				// 0.0 = normal saturation. Higher values mean more saturated image.
  #define GAMMA 1.00f						//1.0 is default brightness. lower values will brighten image, higher values will darken image	

  //#define MOTIONBLUR
	#define MOTIONBLUR_AMOUNT 2.5

  #define WATER_SHADER
	#define WATER_REFLECTIONS
	#define REFLECTION_LENGTH 0.9f
	#define HQ_REFLECTIONS
	#define REFRACT_AMOUNT 0.2f
	#define ABERRATION_AMOUNT 0.105f
	
  #define GLOSSY_REFLECTIONS

  //#define SCREEN_SPACE_RADIOSITY
	#define RADIOSITY_AMOUNT 2.5

#define BANDING_FIX_FACTOR 5.5f


//X-CUSTOM-BEGIN
#define REFLECTIONS
#define QUALITY 10.0
#define RANGE 1.0
//X-CUSTOM-END




// DOF Constants - DO NOT CHANGE
// HYPERFOCAL = (Focal Distance ^ 2)/(Circle of Confusion * F Stop) + Focal Distance
#ifdef USE_DOF
const float HYPERFOCAL = 3.132;
const float PICONSTANT = 3.14159;
#endif





//uniform sampler2D texture;
uniform sampler2D gdepth;
uniform sampler2D gdepthtex;
uniform sampler2D gcolor;
uniform sampler2D gnormal;
uniform sampler2D composite;
uniform sampler2D gaux1;
//uniform sampler2D gaux2;
//uniform sampler2D gaux3; 
//uniform sampler2D gaux4; 

varying vec3 lightVector;

uniform mat4 gbufferProjection;

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

uniform int   isEyeInWater;
uniform float eyeAltitude;
uniform ivec2 eyeBrightness;
uniform ivec2 eyeBrightnessSmooth;

uniform int fogMode;

varying vec4 texcoord;



//Raining
float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;
float wetx  = clamp(wetness, 0.0f, 1.0f);

//Water mask
//float iswater = texture2D(gaux1, texcoord.st).g;
float isice   = 0.0f;

float specularity = texture2D(composite, texcoord.st).g;


// Standard depth function.
float getDepth(vec2 coord) {
    return 2.0 * near * far / (far + near - (2.0 * texture2D(gdepth, coord).x - 1.0) * (far - near));
}
float getWaterDepth(vec2 coord) {
    return 2.0 * near * far / (far + near - (2.0 * texture2D(gdepthtex, coord).x - 1.0) * (far - near));
}
float eDepth(vec2 coord) {
	return texture2D(gdepth, coord).x;
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
	
	
	float getWaterSunOcclusion(vec2 coord) {
		float tex = texture2D(composite, coord).r;								//Call the texture carrying material
			  tex *= 255.0f; 												//Scale material info back into 0-255
			  tex = floor(tex);												//Round materials down to make sure they're integers
		
		float isMaterial;													//Create boolean for material mask
		
		if (tex == 8.0f){
			isMaterial = 1.0f;
		} else {
			isMaterial = 0.0f;
		}
		
		return isMaterial;
	}


//Calculate Time of Day

	float timefract = worldTime;
	float timePow = 3.0f;

	float TimeSunrise  = ((clamp(timefract, 23000.0, 24000.0) - 23000.0) / 1000.0) + (1.0 - (clamp(timefract, 0.0, 6000.0)/6000.0));
		  
	float TimeNoon     = ((clamp(timefract, 0.0, 6000.0)) / 6000.0) - ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0);
	  
	float TimeSunset   = ((clamp(timefract, 6000.0, 12000.0) - 6000.0) / 6000.0) - ((clamp(timefract, 12000.0, 12750.0) - 12000.0) / 750.0);
		  
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
			
				
				float sample = 0.0f;

					sample = 1.0f - texture2D(composite, tx + vec2(noise*delta.x, noise*delta.y)).b;
					sample += 1.0f - texture2D(composite, tx + vec2(noise2*delta.x, noise2*delta.y)).b;
					sample += 1.0f - texture2D(composite, tx + vec2(noise3*delta.x, noise3*delta.y)).b;
					sample += 1.0f - texture2D(composite, tx + vec2(noise4*delta.x, noise4*delta.y)).b;
					sample += 1.0f - texture2D(composite, tx + vec2(noise5*delta.x, noise5*delta.y)).b;
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

#ifdef LENS
	
	vec3 DoLens(float scale, vec3 color, float power, float offset, float curve){
			vec3 sP = sunPosition;
			vec3 c;

			vec2 lPos = sP.xy / -sP.z;
			//lPos.y *= aspectRatio;
			lPos.x *= 1.0f/aspectRatio;
			lPos.xy *= 1.40f;						
			lPos = (lPos + 1.0f)/2.0f;
			
			

			vec2 flare1scale = vec2(1.7f*scale, 1.7f*scale);
			float flare1pow = 12.0f;
			vec2 flare1pos = vec2((1.0f - lPos.x)*(offset + 1.0f) - (offset * 0.5f), (1.0f - lPos.y)*(offset + 1.0f) - (offset * 0.5f)) * vec2(aspectRatio*flare1scale.x, flare1scale.y);
			
			
			//float flare1 = distance(flare1pos, texcoord.st);
			float flare1 = distance(flare1pos, vec2(texcoord.s*aspectRatio*flare1scale.x, texcoord.t*flare1scale.y));
												
			
				  flare1 = 0.5 - flare1;
				  flare1 = clamp(flare1, 0.0, 10.0) * clamp(-sP.z, 0.0, 1.0);
				  //flare1 *= sunmask;
				  flare1 = pow(flare1, curve);
				  
				  flare1 *= flare1pow;
				  
				  	c = flare1 * color * power;
		return c;			
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



//fake albedo
vec3 DoFakeAlbedo(vec3 color) {
	color.rgb *= 90.0f;
	
	//Properties
		float tonemapContrast 		= 0.1f;
		float tonemapSaturation 	= 6.2f; 
		float tonemapDecay			= 400.0f;
		float tonemapCurve			= 0.1f;
		

	//color.xyz*=EBrightnessV2;
	color.rgb += 0.001f;
	
	vec3 colorN = normalize(color.rgb);
	
	vec3 clrfr = color.rgb/colorN.rgb;
	     clrfr = pow(clrfr.rgb, vec3(tonemapContrast));
		 
	colorN.rgb = pow(colorN.rgb, vec3(tonemapSaturation));
	
	color.rgb = clrfr.rgb * colorN.rgb;

	return (color.rgb * (1.0 + color.rgb/tonemapDecay))/(color.rgb + tonemapCurve);
}



//X-CUSTOM-ADDON
#ifdef REFLECTIONS

float fresnelPower = 3.0f;
float subpixelRoughness = getWaterDepth(texcoord.st) / far;
float waterWaviness = 25.0f;

vec3 convertScreenSpaceToWorldSpace(vec2 co) {
    vec4 fragposition = gbufferProjectionInverse * vec4(vec3(co, texture2D(gdepthtex, co).x) * 2.0 - 1.0, 1.0);
    fragposition /= fragposition.w;
    return fragposition.xyz;
}

vec3 convertScreenSpaceToWorldSpaceFake(vec2 co) {
    vec4 fragposition = gbufferProjectionInverse * vec4(vec3(co, 0.0f) * 2.0 - 1.0, 1.0);
    fragposition /= fragposition.w;
    return fragposition.xyz;
}

vec3 convertCameraSpaceToScreenSpace(vec3 cameraSpace) {
    vec4 clipSpace = gbufferProjection * vec4(cameraSpace, 1.0);
    vec3 NDCSpace = clipSpace.xyz / clipSpace.w;
    vec3 screenSpace = 0.5 * NDCSpace + 0.5;
		 screenSpace.z = 0.1f;
    return screenSpace;
}



vec4 ComputeFakeSkyReflection(vec3 col, vec2 waves) {

    vec3 cameraSpacePosition = convertScreenSpaceToWorldSpace(texcoord.st);
    vec3 cameraSpaceNormal = texture2D(gnormal, texcoord.st).rgb * 2.0f - 1.0f;
		 cameraSpaceNormal -= vec3(waves, 0.0f) * waterWaviness;
    vec3 cameraSpaceViewDir = normalize(cameraSpacePosition);
	vec4 color = vec4(0.0f);
	
	float skylight = 1.0f;
	   color.rgb = col * skylight;
	   color.a   = pow(clamp(1.0f + dot(cameraSpaceViewDir, cameraSpaceNormal), 0.0f, 1.0f), fresnelPower) * 0.8f + 0.1f;
	   
	return color;
	
}

vec4 ComputeFakeReflection(float noise1, float noise2, float noise3, vec2 waves) {
	
	float stochasticAmount = 0.05f * subpixelRoughness;
    float initialStepAmount = 1.0 - clamp(0.01 / 100.0, 0.0, 0.99);
		  initialStepAmount *= 50.0f;
	float stepRefinementAmount = 0.7f;
	int maxRefinements = 0;
	
    vec2 screenSpacePosition2D = texcoord.st;
    vec3 cameraSpacePosition = convertScreenSpaceToWorldSpaceFake(screenSpacePosition2D);
	
    vec3 cameraSpaceNormal = texture2D(gnormal, screenSpacePosition2D).rgb * 2.0f - 1.0f;
		 cameraSpaceNormal += vec3(noise1, noise2, noise3) * stochasticAmount;
		 cameraSpaceNormal += cameraSpaceNormal * noise1 * stochasticAmount * 0.0f;
		 cameraSpaceNormal += vec3(waves, 0.0f) * 5.0f;
		 
    vec3 cameraSpaceViewDir = normalize(cameraSpacePosition);
    vec3 cameraSpaceVector = initialStepAmount * normalize(reflect(cameraSpaceViewDir,cameraSpaceNormal));
	vec3 oldPosition = cameraSpacePosition;
    vec3 cameraSpaceVectorPosition = oldPosition + cameraSpaceVector;
    vec3 currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
    vec4 color = vec4(texture2D(gcolor, screenSpacePosition2D).rgb, 0.0);
	int numRefinements = 0;
    int count = 0;

   // while(count < far/initialStepAmount*RANGE)
   // {
       // if(currentPosition.x < 0 || currentPosition.x > 1 ||
       //    currentPosition.y < 0 || currentPosition.y > 1 ||
       //    currentPosition.z < 0 || currentPosition.z > 1) { 

		//   break;
		   
		//   }

        vec2 samplePos = currentPosition.xy;
        float sampleDepth = convertScreenSpaceToWorldSpaceFake(samplePos).z;
        float currentDepth = cameraSpaceVectorPosition.z;
        float diff = sampleDepth - currentDepth;
        float error = length(cameraSpaceVector);
        if(diff >= 0 && diff <= error) {
			cameraSpaceVector *= stepRefinementAmount;
			cameraSpaceVectorPosition = oldPosition;
			numRefinements++;
			if(numRefinements >= maxRefinements){
				vec3 normalAtPos = texture2D(gnormal, samplePos).xyz * 2.0 - 1.0;
				float orientation = dot(cameraSpaceVector,normalAtPos);
				color = texture2D(gcolor, samplePos);
				color.a *= pow(clamp(1 + dot(cameraSpaceViewDir,cameraSpaceNormal), 0.0, 1.0), fresnelPower);
				color.a *= clamp(1 - pow(distance(vec2(0.5), samplePos)*2.0, 2.0), 0.0, 2.0);
				//break;
			}
		}
		
		oldPosition = cameraSpaceVectorPosition;
        cameraSpaceVectorPosition += cameraSpaceVector;
		currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
       // count++;
   // }
	
    return color * vec4(1.0f, 1.0f, 1.0f, 1.0f);
	
}

vec4 ComputeReflection(float noise1, float noise2, float noise3, float specularity) {
	float reflectionRange = RANGE * 0.10f;
	float stochasticAmount = 0.0f;
    float initialStepAmount = 1.0 - clamp(0.01f / 100.0, 0.0, 0.99) + noise1 *  0.0f;
	float stepRefinementAmount = .1;
	int maxRefinements = 0;
	
    vec2 screenSpacePosition2D = texcoord.st;
    vec3 cameraSpacePosition = convertScreenSpaceToWorldSpace(screenSpacePosition2D);
	
    vec3 cameraSpaceNormal = texture2D(gnormal, screenSpacePosition2D).rgb * 2.0f - 1.0f;
		 cameraSpaceNormal += vec3(noise1, noise2, noise3) * stochasticAmount;

    vec3 cameraSpaceViewDir = normalize(cameraSpacePosition);
    vec3 cameraSpaceVector = initialStepAmount * normalize(reflect(cameraSpaceViewDir,cameraSpaceNormal));
	vec3 oldPosition = cameraSpacePosition;
    vec3 cameraSpaceVectorPosition = oldPosition + cameraSpaceVector;
    vec3 currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
    vec4 color = vec4(texture2D(gcolor, screenSpacePosition2D).rgb, 0.0);
	int numRefinements = 0;
    int count = 0;

    while(count < far/initialStepAmount*reflectionRange)
    {
        if(currentPosition.x < 0 || currentPosition.x > 1 ||
           currentPosition.y < 0 || currentPosition.y > 1 ||
           currentPosition.z < 0 || currentPosition.z > 1) { 

		   break;
		   
		   }

        vec2 samplePos = currentPosition.xy;
        float sampleDepth = convertScreenSpaceToWorldSpace(samplePos).z;
        float currentDepth = cameraSpaceVectorPosition.z;
        float diff = sampleDepth - currentDepth;
        float error = length(cameraSpaceVector);
        if(diff >= 0 && diff <= error) {
			cameraSpaceVector *= stepRefinementAmount;
			cameraSpaceVectorPosition = oldPosition;
			numRefinements++;
			//if(numRefinements >= maxRefinements){
				color = texture2D(gcolor, samplePos);
				color.a *= mix(pow(clamp(1 + dot(cameraSpaceViewDir,cameraSpaceNormal), 0.0, 2.0), 3.0f), 1.0f, 0.1f);
				color.a *= clamp(1 - pow(distance(vec2(0.5), samplePos)*2.0, 2.0), 0.0, 1.0);
				break;
			//}
		}
		
		cameraSpaceVector *= 2.5f;	//Each step gets bigger
		
		oldPosition = cameraSpaceVectorPosition;
        cameraSpaceVectorPosition += cameraSpaceVector;
		currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
        count++;
    }
	
	
    return color;
}

vec4 ComputeWaterReflection(float noise1, float noise2, float noise3, vec2 waves, vec2 offset) {
	float reflectionRange = RANGE * 0.45f;
	float stochasticAmount = 0.0f * subpixelRoughness + 0.00f;
    float initialStepAmount = 1.0 - clamp(0.01f / 100.0, 0.0, 0.99) + noise1 *  0.0f;
		  initialStepAmount *= 2.0f;
	float stepRefinementAmount = .1;
	int maxRefinements = 0;
	
    vec2 screenSpacePosition2D = texcoord.st + offset * 1.0f;
    vec3 cameraSpacePosition = convertScreenSpaceToWorldSpace(screenSpacePosition2D);
	
    vec3 cameraSpaceNormal = texture2D(gnormal, screenSpacePosition2D).rgb * 2.0f - 1.0f;
		 cameraSpaceNormal += vec3(noise1, noise2, noise3) * stochasticAmount;
		 cameraSpaceNormal += cameraSpaceNormal * noise1 * stochasticAmount * 0.0f;
		 cameraSpaceNormal -= vec3(waves, -(1.0f - normalize(distance(waves, vec2(0.0f)))) * 0.0f) * waterWaviness;
		 cameraSpaceNormal += vec3(offset, 0.0f) * 0.0f;
		 
    vec3 cameraSpaceViewDir = normalize(cameraSpacePosition);
    vec3 cameraSpaceVector = initialStepAmount * normalize(reflect(cameraSpaceViewDir,cameraSpaceNormal));
	vec3 oldPosition = cameraSpacePosition;
    vec3 cameraSpaceVectorPosition = oldPosition + cameraSpaceVector;
    vec3 currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
    vec4 color = vec4(texture2D(gcolor, screenSpacePosition2D).rgb, 0.0);
	int numRefinements = 0;
    int count = 0;

    while(count < far/initialStepAmount*reflectionRange)
    {
        if(currentPosition.x < 0 || currentPosition.x > 1 ||
           currentPosition.y < 0 || currentPosition.y > 1 ||
           currentPosition.z < 0 || currentPosition.z > 1) { 

		   break;
		   
		   }

        vec2 samplePos = currentPosition.xy;
        float sampleDepth = convertScreenSpaceToWorldSpace(samplePos).z;
        float currentDepth = cameraSpaceVectorPosition.z;
        float diff = sampleDepth - currentDepth;
        float error = length(cameraSpaceVector);
        if(diff >= 0 && diff <= error) {
			cameraSpaceVector *= stepRefinementAmount;
			cameraSpaceVectorPosition = oldPosition;
			numRefinements++;
			//if(numRefinements >= maxRefinements){
				color = texture2D(gcolor, samplePos);
				//color.a *= pow(clamp(1 + dot(cameraSpaceViewDir,cameraSpaceNormal), 0.0, 1.0), fresnelPower);
				color.a *= clamp(1 - pow(distance(vec2(0.5), samplePos)*2.0, 2.0), 0.0, 1.0);
				break;
			//}
		}
		
		cameraSpaceVector *= 1.25f;	//Each step gets bigger
		
		oldPosition = cameraSpaceVectorPosition;
        cameraSpaceVectorPosition += cameraSpaceVector;
		currentPosition = convertCameraSpaceToScreenSpace(cameraSpaceVectorPosition);
        count++;
    }
	
	/*
	if (count > far/initialStepAmount*RANGE) {
	
		//color = ComputeFakeReflection(noise1, noise2, noise3, waves);
		
	}
	*/
    return color;
}
#endif
//X-CUSTON-END


// Main --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void main() {

	vec4 color = texture2D(gcolor, texcoord.st);
	


	//Attempt to hide LDR artifacts
		color.r = pow(color.r, BANDING_FIX_FACTOR);
		color.g = pow(color.g, BANDING_FIX_FACTOR);
		color.b = pow(color.b, BANDING_FIX_FACTOR);
	
	int iswater;
	
	if(texture2D(gaux1, texcoord.st).r > 0.15f){
		iswater = 1;
	} else {
		iswater = 0;
	}
	
	//Land/sky mask
	float land = getLand(texcoord.st);
	
	//Eye sky factor
	float eyeSkylightFactor = eyeBrightnessSmooth.y / 16.0f;
	
		  eyeSkylightFactor = min(eyeSkylightFactor, 16.0f) / 16.0f;
	
	
//Curve times
		  TimeSunrise  = pow(TimeSunrise, timePow);
		  TimeNoon     = pow(TimeNoon, 1.0f/timePow);
		  TimeSunset   = pow(TimeSunset, timePow);
		  TimeMidnight = pow(TimeMidnight, 1.0f/timePow);
	
//Common variables

	float depth = eDepth(texcoord.xy);
	vec2 Texcoord2 = texcoord.st;
	float linDepth = getDepth(texcoord.st);
	vec3 normal = texture2D(gnormal, texcoord.st).rgb;
	vec3 normalBiased = normal * 2.0f - 1.0f;
	
	
	
//Fragposition
	vec4 currentPosition = vec4(texcoord.x * 2.0f - 1.0f, texcoord.y * 2.0f - 1.0f, 2.0f * depth - 1.0f, 1.0f);

	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0f - 1.0f, texcoord.t * 2.0f - 1.0f, 2.0f * texture2D(gdepth, texcoord.st).x - 1.0f, 1.0f);
	fragposition /= fragposition.w;
	
	vec3 npos = normalize(fragposition.xyz);
	
	vec3 specular = reflect(npos, normalBiased);


	
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
	
	vec3 sunrise_amb;
	 sunrise_amb.r = 0.00 * TimeSunrise;
	 sunrise_amb.g = 0.23 * TimeSunrise;
	 sunrise_amb.b = 0.999 * TimeSunrise;	 
	 
	
	vec3 noon_sun;
	 noon_sun.r = mix(1.00, 1.00, rayleigh) * TimeNoon;
	 noon_sun.g = mix(1.00, 0.48, rayleigh) * TimeNoon;
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
	 sunset_sun.g = 0.38 * TimeSunset;
	 sunset_sun.b = 0.0 * TimeSunset;
	
	vec3 sunset_amb;
	 sunset_amb.r = 0.252 * TimeSunset;
	 sunset_amb.g = 0.427 * TimeSunset;
	 sunset_amb.b = 0.999 * TimeSunset;
	
	vec3 midnight_sun;
	 midnight_sun.r = 0.45 * 0.8 * 0.325 * TimeMidnight;
	 midnight_sun.g = 0.6 * 0.8 * 0.325 * TimeMidnight;
	 midnight_sun.b = 0.8 * 0.8 * 0.325 * TimeMidnight;
	
	vec3 midnight_amb;
	 midnight_amb.r = 0.3 * 0.005 * TimeMidnight;
	 midnight_amb.g = 0.4 * 0.005 * TimeMidnight;
	 midnight_amb.b = 0.8 * 0.005 * TimeMidnight;


	vec3 sunlight_color;
	 sunlight_color.r = sunrise_sun.r + noon_sun.r + sunset_sun.r + midnight_sun.r;
	 sunlight_color.g = sunrise_sun.g + noon_sun.g + sunset_sun.g + midnight_sun.g;
	 sunlight_color.b = sunrise_sun.b + noon_sun.b + sunset_sun.b + midnight_sun.b;
	
	vec3 ambient_color;
	 ambient_color.r = sunrise_amb.r + noon_amb.r + sunset_amb.r + midnight_amb.r;
	 ambient_color.g = sunrise_amb.g + noon_amb.g + sunset_amb.g + midnight_amb.g;
	 ambient_color.b = sunrise_amb.b + noon_amb.b + sunset_amb.b + midnight_amb.b;
	 
	vec3 reflected_color;
	 reflected_color = mix(sunlight_color, ambient_color, 0.5f);
	 //reflected_color = mix(vec3(0.64f, 0.73f, 0.34f), reflected_color, 0.5f);
	 //reflected_color = sunlight_color;
	 
	vec3 ambfill_color;
	 ambfill_color = mix(sunlight_color, ambient_color, 0.55f);
	 
	 //ambient_color = mix(ambient_color, vec3(dot(ambient_color, vec3(1.0))), SKY_DESATURATION);
	 
	vec3 skycolor = mix(sunlight_color, vec3(1.0f), 0.5f);
	 
	 float sun_fill = 0.251f * (1.0f - TimeMidnight);
	
	 ambient_color = mix(ambient_color, sunlight_color, sun_fill);
	 vec3 ambient_color_rain = vec3(0.09, 0.09, 0.09) * (1.0f - TimeMidnight * 0.95f); //rain
	 ambient_color = mix(ambient_color, ambient_color_rain, rainx); //rain
	 
	 //Linearize colors
	 
	 ambient_color = pow(ambient_color, vec3(0.454545f));
	 sunlight_color = pow(sunlight_color, vec3(1.454545f));
	 reflected_color = pow(reflected_color, vec3(0.454545f));
	 ambfill_color = pow(ambfill_color, vec3(0.454545f));


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
					
					
						float width11 = 10.0;
						float height11 = 10.0;
						float noiseX11 = ((fract(1.0-Texcoord2.s*(width11/2.0))*0.25)+(fract(Texcoord2.t*(height11/2.0))*0.75))*2.0-1.0;
						float noiseY11 = ((fract(1.0-Texcoord2.s*(width11/2.0))*0.75)+(fract(Texcoord2.t*(height11/2.0))*0.25))*2.0-1.0;

						
							noiseX11 = clamp(fract(sin(dot(Texcoord2 ,vec2(28.9898,61.633))) * 59468.5453),0.0,1.0)*2.0-1.0;
							noiseY11 = clamp(fract(sin(dot(Texcoord2 ,vec2(26.9898,97.233)*2.0)) * 21388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX11 *= (0.10f*noiseamp);
						noiseY11 *= (0.10f*noiseamp);
						
						float width12 = 11.0;
						float height12 = 11.0;
						float noiseX12 = ((fract(1.0-Texcoord2.s*(width12/2.0))*0.25)+(fract(Texcoord2.t*(height12/2.0))*0.75))*2.0-1.0;
						float noiseY12 = ((fract(1.0-Texcoord2.s*(width12/2.0))*0.75)+(fract(Texcoord2.t*(height12/2.0))*0.25))*2.0-1.0;

						
							noiseX12 = clamp(fract(sin(dot(Texcoord2 ,vec2(30.9898,64.633))) * 61468.5453),0.0,1.0)*2.0-1.0;
							noiseY12 = clamp(fract(sin(dot(Texcoord2 ,vec2(34.9898,99.233)*2.0)) * 23388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX12 *= (0.10f*noiseamp);
						noiseY12 *= (0.10f*noiseamp);		

						float width13 = 12.0;
						float height13 = 12.0;
						float noiseX13 = ((fract(1.0-Texcoord2.s*(width13/2.0))*0.25)+(fract(Texcoord2.t*(height13/2.0))*0.75))*2.0-1.0;
						float noiseY13 = ((fract(1.0-Texcoord2.s*(width13/2.0))*0.75)+(fract(Texcoord2.t*(height13/2.0))*0.25))*2.0-1.0;

						
							noiseX13 = clamp(fract(sin(dot(Texcoord2 ,vec2(32.9898,66.633))) * 63468.5453),0.0,1.0)*2.0-1.0;
							noiseY13 = clamp(fract(sin(dot(Texcoord2 ,vec2(36.9898,101.233)*2.0)) * 25388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX13 *= (0.10f*noiseamp);
						noiseY13 *= (0.10f*noiseamp);		

						float width14 = 13.0;
						float height14 = 13.0;
						float noiseX14 = ((fract(1.0-Texcoord2.s*(width14/2.0))*0.25)+(fract(Texcoord2.t*(height14/2.0))*0.75))*2.0-1.0;
						float noiseY14 = ((fract(1.0-Texcoord2.s*(width14/2.0))*0.75)+(fract(Texcoord2.t*(height14/2.0))*0.25))*2.0-1.0;

						
							noiseX14 = clamp(fract(sin(dot(Texcoord2 ,vec2(34.9898,68.633))) * 65468.5453),0.0,1.0)*2.0-1.0;
							noiseY14 = clamp(fract(sin(dot(Texcoord2 ,vec2(38.9898,103.233)*2.0)) * 27388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX14 *= (0.10f*noiseamp);
						noiseY14 *= (0.10f*noiseamp);	

						float width15 = 14.0;
						float height15 = 14.0;
						float noiseX15 = ((fract(1.0-Texcoord2.s*(width15/2.0))*0.25)+(fract(Texcoord2.t*(height15/2.0))*0.75))*2.0-1.0;
						float noiseY15 = ((fract(1.0-Texcoord2.s*(width15/2.0))*0.75)+(fract(Texcoord2.t*(height15/2.0))*0.25))*2.0-1.0;

						
							noiseX15 = clamp(fract(sin(dot(Texcoord2 ,vec2(36.9898,70.633))) * 67468.5453),0.0,1.0)*2.0-1.0;
							noiseY15 = clamp(fract(sin(dot(Texcoord2 ,vec2(40.9898,105.233)*2.0)) * 29388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX15 *= (0.10f*noiseamp);
						noiseY15 *= (0.10f*noiseamp);	

						float width16 = 15.0;
						float height16 = 15.0;
						float noiseX16 = ((fract(1.0-Texcoord2.s*(width16/2.0))*0.25)+(fract(Texcoord2.t*(height16/2.0))*0.75))*2.0-1.0;
						float noiseY16 = ((fract(1.0-Texcoord2.s*(width16/2.0))*0.75)+(fract(Texcoord2.t*(height16/2.0))*0.25))*2.0-1.0;

						
							noiseX16 = clamp(fract(sin(dot(Texcoord2 ,vec2(38.9898,72.633))) * 69468.5453),0.0,1.0)*2.0-1.0;
							noiseY16 = clamp(fract(sin(dot(Texcoord2 ,vec2(42.9898,107.233)*2.0)) * 31388.5453),0.0,1.0)*2.0-1.0;
						
						noiseX16 *= (0.10f*noiseamp);
						noiseY16 *= (0.10f*noiseamp);		

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
/*
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
	
	col += texture2D(gcolor, texcoord.st).rgb;
						  
						  
						dweight =   dofWeight(dofblur, (vec2( 0.0, 0.4)*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0, 0.4)*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( 0.15,0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur).rgb * dweight;
						
						dweight =   dofWeight(dofblur, (vec2( 0.29,0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( -0.37,0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur).rgb * dweight;

						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur).rgb * dweight;
	
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9).rgb * dweight;	
	
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7).rgb * dweight;
	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4).rgb * dweight;
	
						dweight =   dofWeight(dofblur, (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4);
						dweightall += dweight;
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4).rgb * dweight;	

	color.rgb = col/(dweightall + 0.001);	
	
	
	
	#else
	
	//LQ
	vec4 col = vec4(0.0);
	col += texture2D(gcolor, texcoord.st);
	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur);
	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,-0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,-0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.15,0.37 )*aspectcorrect) * dofblur*0.9);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.37,0.15 )*aspectcorrect) * dofblur*0.9);		
	col += texture2D(gcolor, texcoord.st + (vec2( -0.37,-0.15 )*aspectcorrect) * dofblur*0.9);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.15,-0.37 )*aspectcorrect) * dofblur*0.9);	
	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.7);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.7);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.7);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.7);
			 
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,0.29 )*aspectcorrect) * dofblur*0.4);
	col += texture2D(gcolor, texcoord.st + (vec2( 0.4,0.0 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.29,-0.29 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,-0.4 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,0.29 )*aspectcorrect) * dofblur*0.4);
	col += texture2D(gcolor, texcoord.st + (vec2( -0.4,0.0 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(gcolor, texcoord.st + (vec2( -0.29,-0.29 )*aspectcorrect) * dofblur*0.4);	
	col += texture2D(gcolor, texcoord.st + (vec2( 0.0,0.4 )*aspectcorrect) * dofblur*0.4);	

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

		fragposition.xyz += cameraPosition;
	
		vec4 previousPosition = fragposition;
		previousPosition.xyz -= previousCameraPosition;
		previousPosition = gbufferPreviousModelView * previousPosition;
		previousPosition = gbufferPreviousProjection * previousPosition;
		previousPosition /= previousPosition.w;
	
		vec2 velocity = (currentPosition - previousPosition).st * 0.04f * MOTIONBLUR_AMOUNT;
	
		int samples = 0;
		
		int offsetcount = -2;
		
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
			
			coord = coord + (offsetcount*velocity);

			color += texture2D(gcolor, clamp(coord - vec2(noiseX2*velocity.x, noiseX2*velocity.y), vec2(0.001f), vec2(0.999f)));
			color += texture2D(gcolor, clamp(coord - vec2(noiseY2*velocity.x, noiseY2*velocity.y), vec2(0.001f), vec2(0.999f)));
			color += texture2D(gcolor, clamp(coord - vec2(noiseX4*velocity.x, noiseX4*velocity.y), vec2(0.001f), vec2(0.999f)));
			color += texture2D(gcolor, clamp(coord - vec2(noiseY4*velocity.x, noiseY4*velocity.y), vec2(0.001f), vec2(0.999f)));
			samples += 4;
			
			offsetcount += 1;
			
			coord = Texcoord2;
		
		}
			color = (color/1.0)/samples;
		}
		
		

	
#endif


/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////WATER//////////////////////////////////////////////////////////////////////////
#ifdef WATER_SHADER
	
	const float rspread = 0.30f;						//How long reflections are spread across the screen
	
	float rdepth = texture2D(gdepthtex, texcoord.st).x;


float pix_x = 1.0f / viewWidth;
float pix_y = 1.0f / viewHeight;

	rdepth = pow(rdepth, 1.0f);

const float wnormalclamp = 0.05f;
	
//Detect water surface normals

	//Compare change in depth texture over 1 pixel and return an angle
		float wnormal_x1 = texture2D(gdepthtex, texcoord.st + vec2(pix_x, 0.0f)).x - texture2D(gdepthtex, texcoord.st).x;
		      //wnormal_x1 += texture2D(gdepthtex, texcoord.st + vec2(pix_x*2.0f, 0.0f)).x - texture2D(gdepthtex, texcoord.st).x;
			  //wnormal_x1 *= 0.5f;
		
		float wnormal_x2 = texture2D(gdepthtex, texcoord.st).x - texture2D(gdepthtex, texcoord.st + vec2(-pix_x, 0.0f)).x;			
		      //wnormal_x2 += texture2D(gdepthtex, texcoord.st).x - texture2D(gdepthtex, texcoord.st + vec2(-pix_x*2.0f, 0.0f)).x;
			  //wnormal_x2 *= 0.5f;
		float wnormal_x = 0.0f;
		
		if(abs(wnormal_x1) > abs(wnormal_x2)){
			wnormal_x = wnormal_x2;
		} else {
			wnormal_x = wnormal_x1;
		}
		wnormal_x /= 1.0f - rdepth;	

		wnormal_x = clamp(wnormal_x, -wnormalclamp, wnormalclamp);
		
		wnormal_x *= rspread*1.0f;
		

			  
			  
		float wnormal_y1  = texture2D(gdepthtex, texcoord.st + vec2(0.0f, pix_y*1.0f)).x - texture2D(gdepthtex, texcoord.st).x;
		      //wnormal_y1 += texture2D(gdepthtex, texcoord.st + vec2(0.0f, pix_y*2.0f)).x - texture2D(gdepthtex, texcoord.st).x;
			  //wnormal_y1 *= 0.5f;
		float wnormal_y2  = texture2D(gdepthtex, texcoord.st).x - texture2D(gdepthtex, texcoord.st + vec2(0.0f, -pix_y)).x;		
		      //wnormal_y2 += texture2D(gdepthtex, texcoord.st).x - texture2D(gdepthtex, texcoord.st + vec2(0.0f, -pix_y*2.0f)).x;		
			  //wnormal_y2 *= 0.5f;
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
	float waves = texture2D(gaux1, texcoord.st).g;
	float wavesraw = waves;
		  waves -= 0.5f;
		  waves *= 1.0 - rdepth;
		  waves *= 100.0f;

	//Detect angle of waves by comparing 1 pixel difference and resolving discontinuities
	float wavesdeltax1 = texture2D(gaux1, texcoord.st).g - texture2D(gaux1, texcoord.st + vec2(-pix_x, 0.0f)).g;
	float wavesdeltax2 = texture2D(gaux1, texcoord.st + vec2(pix_x, 0.0f)).g - texture2D(gaux1, texcoord.st).g;
	float wavesdeltax;
	
		if(abs(wavesdeltax1) > abs(wavesdeltax2)){
			wavesdeltax = wavesdeltax2;
		} else {
			wavesdeltax = wavesdeltax1;
		}
		
		wavesdeltax = clamp(wavesdeltax, -0.1f, 0.1f);
		
		wavesdeltax *= 1.0f - rdepth;
		wavesdeltax *= 30.0f;
		  
		  
	float wavesdeltay1 = texture2D(gaux1, texcoord.st).g - texture2D(gaux1, texcoord.st + vec2(0.0f, -pix_y)).g;
	float wavesdeltay2 = texture2D(gaux1, texcoord.st + vec2(0.0f, pix_y)).g - texture2D(gaux1, texcoord.st).g;
	float wavesdeltay = 0.0f;
	
		if(abs(wavesdeltay1) > abs(wavesdeltay2)){
			wavesdeltay = wavesdeltay2;
		} else {
			wavesdeltay = wavesdeltay1;
		}
		wavesdeltay *= 1.0f - rdepth;
		wavesdeltay *= 30.0f;
		
		wavesdeltay = clamp(wavesdeltay, -0.1f, 0.1f);
		  


//Calculate distance of objects behind water
float waterDepthDiff = min(abs(getDepth(texcoord.st) - getWaterDepth(texcoord.st)), 3.0f) / 3.0f;
	  waterDepthDiff = pow(waterDepthDiff, 0.5f);
float refractdist = waterDepthDiff * REFRACT_AMOUNT * 10.0f;

//Perform refraction
float refractamount = 500.1154f*0.35f*refractdist;
float refractamount2 = 0.0214f*0.05f*refractdist;
float refractamount3 = 0.214f*0.15f*refractdist;
float waberration = ABERRATION_AMOUNT;

	vec3 refracted = vec3(0.0f);
	vec3 refractedmask = vec3(0.0f);
	float bigWaveRefract = 600.0f * (1.0f - rdepth)*refractdist;
	float bigWaveRefractScale = 1000.0f * (1.0f - rdepth)*refractdist;
	
	vec2 bigRefract = vec2(wnormal_x*bigWaveRefract, wnormal_y*bigWaveRefract);
	
	vec2 refractcoord_r = texcoord.st;
	vec2 refractcoord_g = texcoord.st;
	vec2 refractcoord_b = texcoord.st;
	
	for (int i = 0; i < 1; ++i) {
			
				if(iswater != 1.0f) {
					break;
				}
	
			 refractcoord_r = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y) * (waberration * 2.0f + 1.0f);
			 //refractcoord_g = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y) * (waberration + 1.0f);
			 //refractcoord_b = texcoord.st * (1.0f + waves*refractamount3) - (waves*refractamount3/2.0f) + vec2(wavesdeltax*refractamount*(-wnormal_x*0.3f) + waves*refractamount2 + (-wnormal_x*0.4f) - bigRefract.x, wavesdeltay*refractamount*(-wnormal_y*0.3f) + waves*refractamount2 + (-wnormal_y*0.4f) - bigRefract.y);
				
				refractcoord_r = refractcoord_r * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				//refractcoord_g = refractcoord_g * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				//refractcoord_b = refractcoord_b * vec2(1.0f - abs(wnormal_x) * bigWaveRefractScale, 1.0f - abs(wnormal_y) * bigWaveRefractScale) + vec2(abs(wnormal_x) * bigWaveRefractScale * 0.5f, abs(wnormal_y) * bigWaveRefractScale * 0.5f);
				
			
			refractcoord_r.s = clamp(refractcoord_r.s, 0.001f, 0.999f);
			refractcoord_r.t = clamp(refractcoord_r.t, 0.001f, 0.999f);	
			
			//refractcoord_g.s = clamp(refractcoord_g.s, 0.001f, 0.999f);
			//refractcoord_g.t = clamp(refractcoord_g.t, 0.001f, 0.999f);
			
			//refractcoord_b.s = clamp(refractcoord_b.s, 0.001f, 0.999f);
			//refractcoord_b.t = clamp(refractcoord_b.t, 0.001f, 0.999f);
			
			
			
			if (refractcoord_r.s > 1.0 || refractcoord_r.s < 0.0 || refractcoord_r.t > 1.0 || refractcoord_r.t < 0.0 /*||
				refractcoord_g.s > 1.0 || refractcoord_g.s < 0.0 || refractcoord_g.t > 1.0 || refractcoord_g.t < 0.0 ||
				refractcoord_b.s > 1.0 || refractcoord_b.s < 0.0 || refractcoord_b.t > 1.0 || refractcoord_b.t < 0.0*/) {
					break;
				}
			
			if (refractcoord_r.st * vec2(iswater) == vec2(0.0f)) {
				break;
			}			
			
			/*
			if (refractcoord_g.st * vec2(iswater) == vec2(0.0f)) {
				break;
			}			
			
			if (refractcoord_b.st * vec2(iswater) == vec2(0.0f)) {
				break;
			}
			*/
			
			//refracted.r = texture2D(gcolor, refractcoord_r).r;
			//refracted.g = texture2D(gcolor, refractcoord_g).g;
			//refracted.b = texture2D(gcolor, refractcoord_b).b;
			
			refracted.rgb = texture2D(gcolor, refractcoord_r).rgb;
			
			refracted.rgb = pow(refracted.rgb, vec3(BANDING_FIX_FACTOR));
			
			
			refractedmask.r = texture2D(gaux1, refractcoord_r).r;
			//refractedmask.g = texture2D(gaux3, refractcoord_g).b;
			//refractedmask.b = texture2D(gaux3, refractcoord_b).b;
	
			}
			
	color.rgb = mix(color.rgb, refracted.rgb, vec3(refractedmask.r));
	//color.g = mix(color.g, refracted.g, refractedmask.g);
	//color.b = mix(color.b, refracted.b, refractedmask.b);

	
float depthRefracted = getDepth(mix(texcoord.st, refractcoord_r.st, vec2(refractedmask.r)));
float waterDepthRefracted = getWaterDepth(mix(texcoord.st, refractcoord_r.st, vec2(refractedmask.r)));
	
float waterDepthDiffRefract = min( abs( depthRefracted - waterDepthRefracted ), 20.0f ) / 20.0f;

			
float wfresnel = pow(distance(vec2(wnormal_x, wnormal_y) + vec2(wavesdeltax, wavesdeltay) * 0.1, vec2(0.0f)), 0.7f) * 15.0f;

	  
	  
	  
	  
	  
//REFLECTION

	#ifdef WATER_REFLECTIONS
		//color.rgb = worldposition.g;
		
		vec3 reflection = vec3(0.0f);
		float rtransy = 0.01f * rspread;
		float rtransin = 0.05f;
		
		//Water Fog Properties
			//vec3 waterColor = vec3(0.1f, 0.2f, 0.35f);
			vec3  waterColor 			= vec3(60.0f, 50.0f, 40.0f)/255.0f;
				  waterColor 			= pow(waterColor, vec3(1.0f/BANDING_FIX_FACTOR));
			float waterColorDecay		= 0.175f;
			//float waterLightmap 		= mix(texture2D(gnormal, texcoord.st).g, texture2D(gaux2, texcoord.st).b, 0.5f);
			float waterMaxFogDensity	= 0.00f;
			float waterFogDensity 		= mix(0.9f, 0.0f, isEyeInWater);
			float waterFog				= min((waterDepthDiffRefract * waterFogDensity), waterMaxFogDensity);
			waterColor		   		   -= waterFog * waterColorDecay;

		
		#ifdef HQ_REFLECTIONS
		const float rstrong = 3.9f;
		#else
		const float rstrong = 4.0f;
		#endif
		const float reflectwaviness = 0.00095f;
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
		
		
				#ifdef HQ_REFLECTIONS
				
					for (int i = 0; i < 16; ++i) {
				
						if(iswater != 1.0f) {
							samples += 1;
							break;
						}
						
						rcoord = coordnormal + vec2(dwavesdeltax*4.0f + wnormal_x, dwavesdeltay*9.0f + wnormal_y)*(samples * samples - 1)*0.3f*REFLECTION_LENGTH;
						
						thisdepth = texture2D(gdepth, clamp(texcoord.st + rcoord, 0.001f, 0.999f)).x;
						
						depthcheck = (rdepth - thisdepth);
						depthcheck = 1.0f - depthcheck;
						depthcheck = clamp(depthcheck * 140.0 - 139.0f, 0.0f, 1.0f);
						depthcheck2 = clamp(depthcheck * 70.0 - 69.0f, 0.0f, 1.0f);

							reflectmask   = ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck), 0.001f, 0.999f)).r) * ((17 - samples)/17.0f));
							reflection  += 	((pow(texture2D(gcolor, clamp(texcoord.st + (rcoord*depthcheck), 0.001f, 0.999f)).rgb, vec3(BANDING_FIX_FACTOR)) * reflectmask) * ((17 - samples)/17.0f));

						
						reflectmaskhold += reflectmask;

						samples += 1;
					}
				
				#else
				
					for (int i = 0; i < 8; ++i) {
				
						if(iswater != 1.0f) {
							samples += 1;
							break;
						}
						
						rcoord = coordnormal + vec2(dwavesdeltax*4.0f + wnormal_x, dwavesdeltay*4.0f + wnormal_y)*(samples * samples - 1)*0.3f*REFLECTION_LENGTH;
						
						thisdepth = texture2D(gdepth, clamp(texcoord.st + rcoord, 0.001f, 0.999f)).x;
						
						depthcheck = (rdepth - thisdepth);
						depthcheck = 1.0f - depthcheck;
						depthcheck = clamp(depthcheck * 140.0 - 139.0f, 0.0f, 1.0f);
						depthcheck2 = clamp(depthcheck * 70.0 - 69.0f, 0.0f, 1.0f);

							reflectmask   = ((1.0 - texture2D(gaux1, clamp(texcoord.st + (rcoord*depthcheck), 0.001f, 0.999f)).r) * ((9 - samples)/9.0f));
							reflection  += 	((pow(texture2D(gcolor, clamp(texcoord.st + (rcoord*depthcheck), 0.001f, 0.999f)).rgb, vec3(BANDING_FIX_FACTOR)) * reflectmask) * ((9 - samples)/9.0f));

						
						reflectmaskhold += reflectmask;

						samples += 1;
					}
				
				#endif
				
				reflection /= samples - 1;
				reflectmaskhold /= samples - 1;
				
				reflectmaskhold = pow(reflectmaskhold, 1.0f)*2.5f;
				

				/*
				//compensate for dark water plane
				color.rgb *= mix(1.0f, 1.3f, iswater);
							
				//Darken objects behind water
				color.rgb = mix(color.rgb, vec3(color.r * (1.0f - wfresnel), color.g * (1.0f - wfresnel * 0.99f), color.b * (1.0f - wfresnel * 0.98f)) * (1.0 - reflectmaskhold), iswater);
				
				//Water Fog with conservation of energy considered
				color.rgb = mix(color.rgb, waterColor * (1.0 - reflectmaskhold), iswater * waterFog);
				
				//Add reflections to water only
				reflection *= iswater;
				
				color.rgb = color.rgb + (reflection * rstrong);
				*/

	#endif
	
	
	color.rgb *= mix(1.0f, 0.9f, iswater);
	
	//New Reflections
	if (iswater >= 0.5f) {
	
		vec2 offset1 = vec2(0.5f, 0.0f) / viewWidth;
		vec2 offset2 = vec2(-0.5f, 0.0f) / viewWidth;
		vec2 offset3 = vec2(0.0f, 0.5f) / viewHeight;
		vec2 offset4 = vec2(0.0f, -0.5f) / viewHeight;
	
		vec4 waterReflections  = ComputeWaterReflection(noiseX2, noiseX5, noiseX7, vec2(wavesdeltax, wavesdeltay) * 0.9f, offset1);
			 //waterReflections += ComputeWaterReflection(noiseY2, noiseY5, noiseY7, vec2(wavesdeltax, wavesdeltay) * 0.7f, offset2);
			 //waterReflections += ComputeWaterReflection(noiseX4, noiseX6, noiseX8, vec2(wavesdeltax, wavesdeltay) * 0.5f, offset3);
			 //waterReflections += ComputeWaterReflection(noiseY4, noiseY6, noiseY8, vec2(wavesdeltax, wavesdeltay) * 0.3f, offset4);
			 
			 //waterReflections += ComputeWaterReflection(noiseX9, noiseX11, noiseX13, vec2(wavesdeltax, wavesdeltay) * 0.9f, offset1);
			 //waterReflections += ComputeWaterReflection(noiseY9, noiseY11, noiseY13, vec2(wavesdeltax, wavesdeltay) * 0.7f, offset2);
			 //waterReflections += ComputeWaterReflection(noiseX10, noiseX12, noiseX14, vec2(wavesdeltax, wavesdeltay) * 0.5f, offset3);
			 //waterReflections += ComputeWaterReflection(noiseY10, noiseY12, noiseY14, vec2(wavesdeltax, wavesdeltay) * 0.3f, offset4);
		
		vec4 reflectionsAverage = waterReflections / 1.0f;
		
		float waterSkylight = max(0.0f, texture2D(gaux1, texcoord.st).r * 2.0f - 1.0f);
			  waterSkylight *= 1.0f;
		vec3 fakeSkyColor = pow(ambient_color, vec3(0.2f)) * 0.9f * waterSkylight * (1.0f - isEyeInWater);
		
		vec4 fakeSkyReflection = ComputeFakeSkyReflection(fakeSkyColor, vec2(wavesdeltax, wavesdeltay));
			 
			 reflectionsAverage.rgb = mix(fakeSkyReflection.rgb, reflectionsAverage.rgb, pow(reflectionsAverage.a, 0.5f));
	
		//color.rgb = mix(color.rgb, pow(reflectionsAverage.rgb, vec3(BANDING_FIX_FACTOR)), reflectionsAverage.a);
		color.rgb = mix(color.rgb, pow(reflectionsAverage.rgb, vec3(BANDING_FIX_FACTOR)), fakeSkyReflection.a);
		//color.rgb = vec3(waterSkylight);
	}
		
		float waterRoughness = clamp(1.0f - texture2D(gdepthtex, texcoord.st).x, 0.0f, 1.0f);
		float waterRoughness2 = clamp(texture2D(gdepthtex, texcoord.st).x, 0.0f, 1.0f);
		
		vec3 wnormalsLocal = normalize(normalBiased.rgb + vec3(wavesdeltax, wavesdeltay, 0.0f) * (40.0f * waterRoughness2));
		
		float NdotL = max(dot(wnormalsLocal, lightVector), 0.0f);
		vec3 halfVector;
		float HdotN;
		
			
				halfVector = normalize(lightVector - normalize(fragposition.xyz));
				HdotN = max(0.0f, dot(halfVector, wnormalsLocal));
			
		
		//Water specular highlight
		float waterSpec = HdotN;
		  waterSpec = pow(waterSpec, 111151.0f * waterRoughness + 50.0f);
		  
		  if (waterSpec < 0.0f * pow(waterRoughness, 0.5f)) {
			waterSpec = 0.0f;
		  }
		  
		  waterSpec *= 5.01f * (wfresnel + 0.2);
		  waterSpec *= mix(1.0f, 0.0f, rainx);
		  waterSpec *= mix(0.0f, 1.0f, iswater);
		  waterSpec *= getWaterSunOcclusion(texcoord.st);
		  
		  
		  
		 color.rgb += waterSpec * sunlight_color;
		 
#endif

//color.rgb = wnormals * 0.5f;


#ifdef GLOSSY_REFLECTIONS

//color.rgb = texture2D(composite, texcoord.st).rgb;

/*
const float glosslength = 0.020; //0.095
const float glosslength2 = 0.00575f;
const float g_distance = 50.0f;
const float fadefactor = 0.10f;
float gweight = 0.0f;
float gweight_add = 0.0f;
float compositepass = 0.0f;
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
*/
//float gfresnel = pow(distance(normalBiased, vec3(0.0, 0.0, 1.0))*1.5f, 1.0f);
			
			//vec3 fakeAlbedo = DoFakeAlbedo(color.rgb);
	  
	  		float wetmask = clamp(texture2D(gaux1, texcoord.st).b - 0.25f, 0.0f, 0.1f) / 0.1f * wetx;
			//float distmask = clamp(g_distance * fadefactor - linDepth * fadefactor, 0.0f, 1.0f);
			float g_spec = specularity;
			float g_irr = 0.0f;
			
			float totalspec = g_spec + g_irr;

#ifdef REFLECTIONS
if(totalspec > 0.0f){
	vec4 reflected = ComputeReflection(noiseX2, noiseY2, noiseX4, specularity);
	color.rgb = mix(color.rgb, pow(reflected.rgb, vec3(BANDING_FIX_FACTOR)), reflected.a * 1.0f * totalspec);
}
#else



/*
		
			for (int i = 0; i < 8; ++i) {
			
				if (linDepth > g_distance || totalspec < 0.001f) {
					break;
				}
				
				if (land == 0.0f) {
					break;
				}
						
				glosscoord += ((vec2(normalBiased.x + (noiseX2 * gnoise * normalBiased.x) * aspectRatio, normalBiased.y * 2.0f + (noiseX2 * gnoise * normalBiased.y * 2.0f)))*glosslength);

				
				gdepthdetect = texture2D(gdepth, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).x;

				gdepthcheck = (gdepthr - gdepthdetect);
				gdepthcheck = 1.0f - gdepthcheck;
				gdepthcheck = clamp(gdepthcheck * 480.0 - 479.0f, 0.0f, 1.0f);
				
				
				compositepass = clamp(distance(texture2D(gnormal, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).rgb, normal.rgb) * 2.5f - 0.5f, 0.0f, 2.0f);
				
				gweight = 1.0f * compositepass * gfresnel * gdepthcheck * (9 - i);
				
				gweight_add += gweight;


				
				gloss += pow(max((texture2D(gcolor, clamp(gcoord.st + glosscoord, 0.001f, 0.999f)).rgb - 0.0f) * gweight, 0.0f), vec3(BANDING_FIX_FACTOR));

			}		
		//}

			float finalweight =  (((gweight_add/105.0)*land) * (1.0f - iswater)) * 1.1;

			gloss /= gweight_add * 8811151.0f;
			gloss = max(gloss, 0.0f);
			
			//reflect
			color.rgb = mix(color.rgb, gloss, clamp(finalweight * distmask * g_spec, 0.0f, 1.0f) * pow(dot(gloss, vec3(1.0f)), 0.2f));

			//metallic reflect
			color.rgb += gloss * albedo.rgb * finalweight * distmask * g_irr * pow(dot(gloss, vec3(1.0f)), 0.2f);

	*/		
#endif
#endif



#ifdef SCREEN_SPACE_RADIOSITY

//color.rgb = texture2D(composite, texcoord.st).rgb;

//vec3 normal = texture2D(composite, texcoord.st).rgb * 2.0f - 1.0f;

const float radlength_r = 0.060;
const float radlength_r2 = 0.00575f;
float gweight_r = 0.0f;
float compositepass_r = 0.0f;
float gdepthcheck_r = 0.0f;
const float gnoise_r = 0.000f;
float gnoise_r2 = 0.125f;
float gdepthdetect_r = 0.0f;
float gweight_racc = 0.0f;


float sphere_2 = 0.0f;

vec3 rad = vec3(0.0f);
vec2 radcoord = vec2(0.0f);
vec2 gcoordr = texcoord.st;

//float gfresnel = pow(distance(normalBiased, vec3(0.0, 0.0, 1.0))*1.5f, 2.0f);
	 // gfresnel = 4.0f;

			for (int i = 0; i < 5; ++i) {
			
				for (int k = 0; k < 6; ++k) {
				
					float rradius = pow(i, 1.5f);
			
					sphere_2 = (-1.0f + k + (mod(i, 2) - 1.0f))*rradius/2.0f;
						
					radcoord = ((vec2(normalBiased.x * rradius + normalBiased.y * (sphere_2) + (noiseX2 * gnoise_r) * aspectRatio, normalBiased.y * rradius * 2.0f + normalBiased.x * (sphere_2) * 2.0f + (noiseY2 * gnoise_r) * 2.0f))*radlength_r) / clamp(linDepth, 0.0f, 10.0f);
				
					gdepthcheck_r = distance(linDepth, getDepth(texcoord.st + radcoord)) / i;
					gdepthcheck_r = 1.0f - gdepthcheck_r;
					gdepthcheck_r = pow(clamp(gdepthcheck_r * 1.0 - 0.0f, 0.0f, 1.0f), 0.5f);
				
				
					compositepass_r = clamp(distance(texture2D(gnormal, gcoordr.st + radcoord).rg, normal.rg) * 2.5f - 0.5f, 0.0f, 2.0f);
				
					gweight_r = 1.0f * compositepass_r * 4.0f * gdepthcheck_r  * pow((7 - i), 1.0f);
				
					gweight_racc += gweight_r;
					
					rad += max(pow((texture2D(gcolor, gcoordr.st + radcoord).rgb - 0.0f), vec3(BANDING_FIX_FACTOR)) * gweight_r, vec3(0.0f));
				
				}
				
				if (land == 0.0f) {
					break;
				}
				

				


			}		


		
			
			float finalweight_r =  (((gweight_racc/25.0)*land) * (1.0f - iswater)) * 0.17 * RADIOSITY_AMOUNT;
			
			
			rad /= gweight_racc;
			rad = max(rad, 0.0f);
	
			rad *= max(1.0f - dot(color.rgb, vec3(1.0f)), 0.1f);
	
			//color.rgb += rad * finalweight_r * albedo * 0.15f;
			color.rgb *= 1.0 + ((rad - 0.000f) * finalweight_r * 1.25f);
			//color.rgb *= 1.0f - finalweight_r;
			
			//color.rgb *= max(1.0f - (finalweight_r * (1.0f - min(rad * 6.0f, 1.0f))), 0.0f);

#endif





#ifdef GODRAYS

	float GR = addGodRays(0.0f, Texcoord2, noiseX3, noiseX4, noiseY4, noiseX2, noiseY2, noiseX5, noiseY5, noiseX6, noiseY6)/2.0;

	float GRr = 1.0 - texture2D(composite, texcoord.st).b;
	
	//GR = mix(GR, 0.0f, rainx);
	
	/*
	float GRs  = 1.0 - texture2D(gaux1, vec2(0.55, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.55, 0.45)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.55)).g;
		  GRs += 1.0 - texture2D(gaux1, vec2(0.45, 0.45)).g;

		  GRs /= 3.0;
	*/


	
	

	
	
	GR = pow(GR, 1.0f)*2.5f;
	
	color.r += pow(GR*sunlight_color.r, 1.0f);
	color.g += pow(GR*sunlight_color.g, 1.0f);
	color.b += pow(GR*sunlight_color.b, 1.0f);
	
	
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
	
	//clr += texture2D(gcolor, texcoord.st);
	
	for (int i = 0; i < 1; ++i) {
	//horizontal (70 taps)

	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(5.0f+bloomnoise,5.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(4.0f+bloomnoise,4.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(3.0f+bloomnoise,3.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(2.0f+bloomnoise,2.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(1.0f+bloomnoise,1.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
		//clr += texture2D(gcolor, texcoord.st + (vec2(0.0f,0.0f))*radius)*6.0f;
		
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-1.0f+bloomnoise,1.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-2.0f+bloomnoise,2.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-3.0f+bloomnoise,3.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-4.0f+bloomnoise,4.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-5.0f+bloomnoise,5.0+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;

	//vertical

	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(5.0+bloomnoise,-5.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(4.0+bloomnoise,-4.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(3.0+bloomnoise,-3.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(2.0+bloomnoise,-2.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(1.0+bloomnoise,-1.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-5.0+bloomnoise,-5.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*1.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-4.0+bloomnoise,-4.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*2.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-3.0+bloomnoise,-3.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*3.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-2.0+bloomnoise,-2.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*4.0f;
	clr +=  max(texture2D(gcolor, texcoord.st + (vec2(-1.0+bloomnoise,-1.0f+bloomnoisey))*radius)*glarex - glaresub, 0.0f)*5.0f;
	
	radius *= 3.0;
	
	clr = (clr/10.0f)/2.0f;
	}
	
	const float clrboost = 0.0;
	
	clr.r = clr.r + (clr.r*(clrboost*2.0)) - (clr.g * clrboost) - (clr.b * clrboost);
	clr.g = clr.g + (clr.g*(clrboost*2.0)) - (clr.r * clrboost) - (clr.b * clrboost);
	clr.b = clr.b + (clr.b*(clrboost*2.0)) - (clr.r * clrboost) - (clr.g * clrboost);

	clr.rgb = pow(clr.rgb, vec3(BANDING_FIX_FACTOR));
	
	color.r = color.r + (clr.r*1.0f)*bloomintensity;
	color.g = color.g + (clr.g*1.0f)*bloomintensity;
	color.b = color.b + (clr.b*1.0f)*bloomintensity;
	color = max(color, 0.0f);
	

#endif




#ifdef GLARE2

	color = color * 1.0f;
	
	float radius2 = 0.006f*GLARE_RANGE;
	const float radius2v = 0.002f;
	const float bloomintensity2 = 0.08f*GLARE_AMOUNT;
	
	const float glarex2 = 2.0f;
	const float glaresub2 = 0.0f;
	
	float bloomnoise2 = noiseX4*0.0f;	
	float bloomnoise2y = noiseY4*0.0f;	

	vec4 clr2 = vec4(0.0f);
	
	//clr2 += texture2D(gcolor, texcoord.st);
	
	//horizontal (70 taps)
	
	for (int i = 0; i < 1; ++i) {

	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(5.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(4.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(3.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(2.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(1.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
		//clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0f,0.0f))*radius2)*6.0f;
		
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(-1.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(-2.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(-3.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(-4.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(-5.0f+bloomnoise2,0.0+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;

	//vertical

	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,-5.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,-4.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,-3.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,-2.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,-1.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,5.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*1.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,4.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*2.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,3.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*3.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,2.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*4.0f;
	clr2 += max(texture2D(gcolor, texcoord.st + (vec2(0.0+bloomnoise2,1.0f+bloomnoise2y))*radius2)*glarex2 - glaresub2, 0.0f)*5.0f;
	
	radius2 /= 3.0;
	
	clr2 = (clr2/10.0f)/2.0f;
	}
	
	const float clr2boost = 0.0;
	
	clr2.r = clr2.r + (clr2.r*(clr2boost*2.0)) - (clr2.g * clr2boost) - (clr2.b * clr2boost);
	clr2.g = clr2.g + (clr2.g*(clr2boost*2.0)) - (clr2.r * clr2boost) - (clr2.b * clr2boost);
	clr2.b = clr2.b + (clr2.b*(clr2boost*2.0)) - (clr2.r * clr2boost) - (clr2.g * clr2boost);

	clr2.rgb = pow(clr2.rgb, vec3(BANDING_FIX_FACTOR));
	
	color.r = color.r + (clr2.r*1.0f)*bloomintensity2;
	color.g = color.g + (clr2.g*1.0f)*bloomintensity2;
	color.b = color.b + (clr2.b*1.0f)*bloomintensity2;
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
			lPos.x *= 1.0f/aspectRatio;
			lPos.xy *= 1.40f;						
			lPos = (lPos + 1.0f)/2.0f;
			//lPos = clamp(lPos, vec2(0.001f), vec2(0.999f));
			

			
			float sunmask = 0.0f;
					

					sunmask += 1.0f - getLand(lPos);
					
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
				
				
				
			//Do flares. DoLens(float scale, vec3 color, float power, float offset, float curve);
				
			color.rgb += DoLens(0.2f, sunlight_color, 29.0f, -2.0f, 7.0f) * sunmask;
			

			color.rgb = clamp(color.rgb, 0.0, 10.0);

#endif



if (fogMode != 0) {
	#ifdef CREPUSCULAR_RAYS
		float crepRays = 0.0f;
	#else
		float crepRays = 0.0f;
	#endif
		
	//Sky Fog
		float skyFogDensity = mix(0.0025f * FOG_DENSITY, 0.010f * FOG_DENSITY, rainx) * land;
			  skyFogDensity *= mix(0.0f, 1.0f, (min(eyeSkylightFactor, 0.25f)/0.25f));
		
		float skyFogFactor  = exp(linDepth * skyFogDensity) - 1.0f;
		vec3  skyFogColor   = pow(ambient_color, vec3(1.0f)) * 1.0f;
			  skyFogColor	*= pow(eyeSkylightFactor, 3.0f);
		vec3  caveFogColor  = vec3(0.01f);
		vec3  fogColor1 	= mix(caveFogColor, skyFogColor, vec3(max(eyeSkylightFactor, 0.75f)/0.25f));
		
		color.rgb = mix(color.rgb, skyFogColor, vec3(clamp(skyFogFactor, 0.0f, 1.0f)));
		
	/*	
	//Haze ambient
		float hazeFogDensity    = 0.0095f * FOG_DENSITY * 0.0f;
		float hazeFogAltitude   = 100.0f - eyeAltitude;
		float hazeFogMaxDensity = 0.021f * FOG_DENSITY;
		
		      hazeFogDensity    = clamp(((-fragposition.g + hazeFogAltitude) * 0.000025f), 0.0f, hazeFogMaxDensity);
		float hazeFogFactor     = exp(linDepth * hazeFogDensity) - 1.0f;
		vec3  hazeFogColor		= reflected_color * 0.0f + crepRays * sunlight_color * 4.0f;;
		
			//color.rgb			= mix(color.rgb, hazeFogColor, vec3(min(hazeFogFactor, 1.0f)));
	*/		  
}




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



	color.rgb *= 3.1f;
	 float TonemapOversat = 54.4f;
	 float TonemapCurve   = 1.57f;
	 
	#ifdef TONEMAP_FILMIC

		TonemapOversat = 115.0f;
		TonemapCurve   = 0.499f;
		
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
	
	#ifdef TONEMAP_FILMIC
	color.rgb *= 1.1f;
	color.rgb -= 0.00f;
	#else
	color.rgb *= 1.5f;
	#endif

#endif
/*
float exposureb = 3.0f;

color.rgb = (1.0f - exp( -color.rgb * exposureb ));
*/

color.rgb *= 1.25f;

#ifdef TONEMAP2
	
	//Pseudo hdr
	eyeSkylightFactor = pow(eyeSkylightFactor, 5.0f);
	eyeSkylightFactor *= mix(1.0f, NIGHT_EXPOSURE_BIAS, TimeMidnight);
	color.rgb /= eyeSkylightFactor * 1.275f + 0.15f;
	
	color.rgb *= 80.0f;
	
	//Properties
		float tonemapContrast 		= 1.1f;
		float tonemapSaturation 	= 1.1f; 
		float tonemapDecay			= 300.0f;
		float tonemapCurve			= 25.0f;
		

	//color.xyz*=EBrightnessV2;
	color.rgb += 0.001f;
	
	vec3 colorN = normalize(color.rgb);
	
	vec3 clrfr = color.rgb/colorN.rgb;
	     clrfr = pow(clrfr.rgb, vec3(tonemapContrast));
		 
	colorN.rgb = pow(colorN.rgb, vec3(tonemapSaturation));
	
	color.rgb = clrfr.rgb * colorN.rgb;

	color.rgb = (color.rgb * (1.0 + color.rgb/tonemapDecay))/(color.rgb + tonemapCurve);

#endif

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
