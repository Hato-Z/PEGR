#version 120

uniform sampler2D gcolor;
uniform sampler2D gdepth;
uniform sampler2D shadow;

//to increase shadow draw distance, edit SHADOWDISTANCE and SHADOWHPL below. Both should be equal. Needs decimal point.
//disabling is done by adding "//" to the beginning of a line.


//ADJUSTABLE VARIABLES

#define HQ              //high quality. Only enable HQ or LQ, not both
//#define LQ            //low quality. Only enable HQ or LQ, not both

#define BLURFACTOR 0.9

#define SHADOWDISTANCE 45.0 

/* SHADOWRES:2048 */
/* SHADOWHPL:45.0 */

//END OF ADJUSTABLE VARIABLES


varying vec4 texcoord;

uniform int worldTime;

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferModelViewInverse;
uniform mat4 shadowProjection;
uniform mat4 shadowModelView;

void main() {
	vec4 fragposition = gbufferProjectionInverse * vec4(texcoord.s * 2.0 - 1.0, texcoord.t * 2.0 - 1.0, 2.0 * texture2D(gdepth, texcoord.st).x - 1.0, 1.0);
	fragposition /= fragposition.w;
	
	#ifdef SHADOWDISTANCE
	float drawdistance = SHADOWDISTANCE;
	float drawdistancesquared = pow(drawdistance, 2);
	#endif

	float distance = sqrt(fragposition.x * fragposition.x + fragposition.y * fragposition.y + fragposition.z * fragposition.z);

	float shading = 1.0;
//
	if (distance < drawdistance && distance > 0.1) {
		// shadows
		vec4 worldposition = gbufferModelViewInverse * fragposition;

		float xzDistanceSquared = worldposition.x * worldposition.x + worldposition.z * worldposition.z;
		float yDistanceSquared  = worldposition.y * worldposition.y;
		
		if (yDistanceSquared < drawdistancesquared) {
			worldposition = shadowModelView * worldposition;
			float comparedepth = -worldposition.z;
			worldposition = shadowProjection * worldposition;
			worldposition /= worldposition.w;
			
			worldposition.st = worldposition.st * 0.5 + 0.5;
				
			if (comparedepth > 0.0 && worldposition.s < 1.0 && worldposition.s > 0.0 && worldposition.t < 1.0 && worldposition.t > 0.0){
				float shadowMult = min(1.0 - xzDistanceSquared / drawdistancesquared, 1.0) * min(1.0 - yDistanceSquared / drawdistancesquared, 1.0);
				float sampleDistance = 1.0 / 2048.0;
					
					
				#ifdef HQ
				
					float ntime = worldTime/24000.0;
					float dawn = clamp(ntime*20.0, 0.0, 1.0) * clamp((1.0-ntime)*20.0, 0.0, 1.0);
					//float dawn -= ;
					//float dusk = clamp((ntime-0.5)*10.0, 0.0, 1.0);
					//float dusk += clamp((0.5-ntime)*10.0-9.0, 0.0, 1.0);
					
					float zoffset = 0.0;
					float offsetx = 0.0003*BLURFACTOR;
					float offsety = 0.0006*BLURFACTOR;
					
					float shadowdarkness = 0.65*dawn;
					float diffthresh = 0.5;
					float bluramount = 0.000055*BLURFACTOR;
					
				
					//float col = 0.0;
					
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, 4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, 3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, -2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-3*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-2*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(3*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(4*bluramount+offsetx, -3*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(5*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(5*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(5*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(5*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 5*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 5*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 5*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 5*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-4*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-4*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-4*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-4*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, -4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, -4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);

					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, -4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/80 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, -4*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);

					//shading = clamp((shading/80 + 0.93), 0.3, 0.7)/0.4 - 1;
					shading = shading/80 + 0.97;
					
				#endif
					
				#ifdef LQ
				
					float zoffset = 0.02;
					float offsetx = 0.0003*BLURFACTOR;
					float offsety = 0.0004*BLURFACTOR;
				
					
					//shadow filtering
					
					float shadowdarkness = 0.65;
					float diffthresh = 0.3;
					float bluramount = 0.00009*BLURFACTOR;
					
					
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 2*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, 0*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(-1*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(0*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(1*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);
					shading += 1.0/16 - shadowMult * (clamp(comparedepth - (0.05 + (texture2D(shadow, worldposition.st + vec2(2*bluramount+offsetx, -1*bluramount+offsety)).z) * (256.0 - 0.05)) - zoffset, 0.0, diffthresh)/diffthresh * shadowdarkness - zoffset);

					shading = shading/16 + 0.86;
					
				#endif
			
			}
		}
	}

	gl_FragData[0] = texture2D(gcolor, texcoord.st);
	gl_FragData[1] = texture2D(gdepth, texcoord.st);
	gl_FragData[3] = vec4(texture2D(gcolor, texcoord.st).rgb * shading, 1.0);
}
