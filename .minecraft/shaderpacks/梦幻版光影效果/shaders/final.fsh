#version 120

/*

Settings by Sonic Ether
More realistic depth-of-field by Azraeil.
God Rays by Blizzard
Bloom shader by CosmicSpore (Modified from original source: http://myheroics.wordpress.com/2008/09/04/glsl-bloom-shader/)
Cross-Processing by Sonic Ether.
High Desaturation effect by Sonic Ether
Shaders 2.0 port of Yourself's Cell Shader, port by an anonymous user.
Bug Fixes by Kool_Kat.

*/

// Place two leading Slashes in front of the following '#define' lines in order to disable an option.
// MOTIONBLUR, HDR, and BOKEH_DOF are very beta shaders. Use at risk of weird results.
// MOTIONBLUR and BOKEH_DOF are not compatable with eachother. Shaders break when you enable both.
// GLARE is still a work in progress, swap it for BLOOM if you don't like the effect.

//#define HDR
//#define BOKEH_DOF
//#define GODRAYS
//#define GODRAYS_EXPOSURE 0.2
//#define GODRAYS_SAMPLES 16
//#define GODRAYS_DECAY 0.9
//#define GODRAYS_DENSITY 0.5
#define GLARE
#define GLARE_AMOUNT 0.55
#define GLARE_RANGE 5
//#define BLOOM
//#define BLOOM_AMOUNT 1.0
//#define BLOOM_RANGE 3
//#define CEL_SHADING
//#define CEL_SHADING_THRESHOLD 0.4
//#define CEL_SHADING_THICKNESS 0.004
#define USE_HIGH_QUALITY_BLUR
#define CROSSPROCESS
#define BRIGHTMULT 1.05                  // 1.0 = default brightness. Higher values mean brighter. 0 would be black.
#define COLOR_BOOST	0.1					// 0.0 = normal saturation. Higher values mean more saturated image.
#define SSAO
#define SSAO_LUMINANCE 0.4
#define SSAO_STRENGTH 1.8               //too much strength causes white highlights on extruding edges and behind objects
#define SSAO_LOOP 6
#define SSAO_MAX_DEPTH 0.5
#define SSAO_SAMPLE_DELTA 0.35
#define MOTIONBLUR
#define HIGHDESATURATE

// DOF Constants - DO NOT CHANGE
// HYPERFOCAL = (Focal Distance ^ 2)/(Circle of Confusion * F Stop) + Focal Distance
#ifdef USE_DOF
const float HYPERFOCAL = 3.132;
const float PICONSTANT = 3.14159;
#endif


//uniform sampler2D texture;
uniform sampler2D gcolor;
uniform sampler2D gdepth;
uniform sampler2D composite;
uniform sampler2D gnormal;
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

varying vec4 texcoord;

// Standard depth function.
float getDepth(vec2 coord) {
    return 2.0 * near * far / (far + near - (2.0 * texture2D(gdepth, coord).x - 1.0) * (far - near));
}
float eDepth(vec2 coord) {
	return texture2D(gdepth, coord).x;
}
float realcolor(vec2 coord) {
	return texture2D(gcolor, coord).r;
}


#ifdef BOKEH_DOF

const float blurclamp = 10.0;  // max blur amount
const float bias = 0.3;	//aperture - bigger values for shallower depth of field

#endif

#ifdef GODRAYS
	vec4 addGodRays(vec4 nc, vec2 tx) {
		float threshold = 0.99 * far;
//		bool foreground = false;
		float depthGD = getDepth(tx);
		if ( (worldTime < 14000 || worldTime > 22000) && (sunPosition.z < 0) && (depthGD < threshold) ) {
			vec2 lightPos = sunPosition.xy / -sunPosition.z;
			lightPos.y *= aspectRatio;
			lightPos = (lightPos + 1.0)/2.0;
			//vec2 coord = tx;
			vec2 delta = (tx - lightPos) * GODRAYS_DENSITY / float(GODRAYS_SAMPLES);
			float decay = -sunPosition.z / 100.0;
			vec3 colorGD = vec3(0.0);
			
			for (int i = 0; i < GODRAYS_SAMPLES; i++) {
				tx -= delta;
				if (tx.x < 0.0 || tx.x > 1.0) {
					if (tx.y < 0.0 || tx.y > 1.0) {
						break;
					}
				}
				vec3 sample = vec3(0.0);
				if (getDepth(tx) > threshold) {
					sample = texture2D(composite, tx).rgb;
				}
				sample *= vec3(decay);
				if (distance(tx, lightPos) > 0.05) {
					sample *= 0.2;
				}
					colorGD += sample;
					decay *= GODRAYS_DECAY;
			}
			return (nc + GODRAYS_EXPOSURE * vec4(colorGD, 0.0));
        } else {
			return nc;
		}
	}
#endif 

/*
#ifdef BLOOM
	vec4 addBloom(vec4 c, vec2 t) {
		int j;
		int i;
		vec4 bloom = vec4(0.0);
		vec2 loc = vec2(0.0);
		float count = 0.0;
		
		for( i= -BLOOM_RANGE ; i < BLOOM_RANGE; i++ ) {
			for ( j = -BLOOM_RANGE; j < BLOOM_RANGE; j++ ) {
				loc = t + vec2(j, i)*0.004;
				
				// Only add to bloom texture if loc is on-screen.
				if(loc.x > 0 && loc.x < 1 && loc.y > 0 && loc.y < 1) {
					bloom += texture2D(composite, loc) * BLOOM_AMOUNT;
					count += 1;
				}
			}
		}
		bloom /= vec4(count);
		
		if (c.r < 0.3)
		{
			return bloom*bloom*0.012;
		}
		else
		{
			if (c.r < 0.5)
			{
				return bloom*bloom*0.009;
			}
			else
			{
				return bloom*bloom*0.0075;
			}
		}
	}
#endif
*/

#ifdef CEL_SHADING
	float getCellShaderFactor(vec2 coord) {
    float d = getDepth(coord);
    vec3 n = normalize(vec3(getDepth(coord+vec2(CEL_SHADING_THICKNESS,0.0))-d,getDepth(coord+vec2(0.0,CEL_SHADING_THICKNESS))-d , CEL_SHADING_THRESHOLD));
    //clamp(n.z*3.0,0.0,1.0);
    return n.z; 
	}
#endif

#ifdef SSAO
uniform float viewWidth;
uniform float viewHeight;

// Alternate projected depth (used by SSAO, probably AA too)
float getProDepth( vec2 coord ) {
	float depth = texture2D(gdepth, coord).x;
	return ( 2.0 * near ) / ( far + near - depth * ( far - near ) );
}

float znear = near; //Z-near
float zfar = far; //Z-far

float diffarea = 0.5; //self-shadowing reduction
float gdisplace = 0.3; //gauss bell center

bool noise = false; //use noise instead of pattern for sample dithering?
bool onlyAO = false; //use only ambient occlusion pass?

vec2 texCoord = texcoord.st;

vec2 rand(vec2 coord) { //generating noise/pattern texture for dithering
  float width = 1.0;
  float height = 1.0;
  float noiseX = ((fract(1.0-coord.s*(width/2.0))*0.25)+(fract(coord.t*(height/2.0))*0.75))*2.0-1.0;
  float noiseY = ((fract(1.0-coord.s*(width/2.0))*0.75)+(fract(coord.t*(height/2.0))*0.25))*2.0-1.0;

  if (noise) {
    noiseX = clamp(fract(sin(dot(coord ,vec2(12.9898,78.233))) * 43758.5453),0.0,1.0)*2.0-1.0;
    noiseY = clamp(fract(sin(dot(coord ,vec2(12.9898,78.233)*2.0)) * 43758.5453),0.0,1.0)*2.0-1.0;
  }
  return vec2(noiseX,noiseY)*0.001;
}

float compareDepths(float depth1, float depth2, int zfar) {  
  float garea = 1.5; //gauss bell width    
  float diff = (depth1 - depth2) * 100.0; //depth difference (0-100)
  //reduce left bell width to avoid self-shadowing 
  if (diff < gdisplace) {
    garea = diffarea;
  } else {
    zfar = 1;
  }

  float gauss = pow(2.7182,-2.0*(diff-gdisplace)*(diff-gdisplace)/(garea*garea));
  return gauss;
} 

float calAO(float depth, float dw, float dh) {  
  float temp = 0;
  float temp2 = 0;
  float coordw = texCoord.x + dw/depth;
  float coordh = texCoord.y + dh/depth;
  float coordw2 = texCoord.x - dw/depth;
  float coordh2 = texCoord.y - dh/depth;

  if (coordw  < 1.0 && coordw  > 0.0 && coordh < 1.0 && coordh  > 0.0){
    vec2 coord = vec2(coordw , coordh);
    vec2 coord2 = vec2(coordw2, coordh2);
    int zfar = 0;
    temp = compareDepths(depth, getProDepth(coord),zfar);

    //DEPTH EXTRAPOLATION:
    if (zfar > 0){
      temp2 = compareDepths(getProDepth(coord2),depth,zfar);
      temp += (1.0-temp)*temp2; 
    }
  }

  return temp;  
}  

vec3 getSSAOFactor() {
	vec2 noise = rand(texCoord); 
	float depth = getProDepth(texCoord);
  if (depth > SSAO_MAX_DEPTH) {
    return vec3(1.0,1.0,1.0);
  }
  float cdepth = texture2D(gdepth,texCoord).g;
	
	float ao;
	float s;
	
  float incx = 1.0 / viewWidth * SSAO_SAMPLE_DELTA;
  float incy = 1.0 / viewHeight * SSAO_SAMPLE_DELTA;
  float pw = incx;
  float ph = incy;
  float aoMult = SSAO_STRENGTH;
  int aaLoop = SSAO_LOOP;
  float aaDiff = (1.0 + 2.0 / aaLoop);
  for (int i = 0; i < aaLoop ; i++) {
    float npw  = (pw + 0.1 * noise.x) / cdepth;
    float nph  = (ph + 0.1 * noise.y) / cdepth;

    ao += calAO(depth, pw, ph) * aoMult;
    ao += calAO(depth, pw, -ph) * aoMult;
    ao += calAO(depth, -pw, ph) * aoMult;
    ao += calAO(depth, -pw, -ph) * aoMult;
    pw += incx;
    ph += incy;
    aoMult /= aaDiff; 
    s += 4.0;
  }
	
	ao /= s;
	ao = 1.0-ao;	
  ao = clamp(ao, 0.0, 0.5) * 2.0;
	
  return vec3(ao);
}
#endif




// Main ---------------------------------------------------------------------------------------------------
void main() {

	vec4 color = texture2D(composite, texcoord.st);


#ifdef BOKEH_DOF

	float depth = eDepth(texcoord.xy);
	
	if (depth > 0.9999) {
		depth = 1.0;
	}
	

	float cursorDepth = eDepth(vec2(0.5, 0.5));
	
	if (cursorDepth > 0.9999) {
		cursorDepth = 1.0;
	}
	
	
	vec2 aspectcorrect = vec2(1.0, aspectRatio) * 1.5;
	
	float factor = (depth - cursorDepth);
	 
	vec2 dofblur = (vec2 (clamp( factor * bias, -blurclamp, blurclamp )))*0.4;
	
	

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

#ifdef GODRAYS
	color.r = addGodRays(color, texcoord.st).r;
	color.g = addGodRays(color, texcoord.st).g;
	color.b = addGodRays(color, texcoord.st).b;
#endif

	


#ifdef MOTIONBLUR

	vec4 aux1   = texture2D(gaux1, texcoord.st);
	vec4 depth  = texture2D(gdepth, texcoord.st);
	

	//gl_FragData[4] = aux1;

	/*
		if (aux1.r > 0.5 || depth.x > 0.9999999) {
		color = texture2D(composite, texcoord.st);
		return;
		}
	*/
	
	/*
	
		if (depth.x < 1.9999999) {
		vec4 currentPosition = vec4(texcoord.x * 2.0 - 1.0, texcoord.y * 2.0 - 1.0, 2.0 * depth.x - 1.0, 1.0);
	
		vec4 fragposition = gbufferProjectionInverse * currentPosition;
		fragposition = gbufferModelViewInverse * fragposition;
		fragposition /= fragposition.w;
		fragposition.xyz += cameraPosition;
	
		vec4 previousPosition = fragposition;
		previousPosition.xyz -= previousCameraPosition;
		previousPosition = gbufferPreviousModelView * previousPosition;
		previousPosition = gbufferPreviousProjection * previousPosition;
		previousPosition /= previousPosition.w;
	
		vec2 velocity = (currentPosition - previousPosition).st * 0.006;
	
		int samples = 1;

	vec2 coord = texcoord.st + velocity;
	for (int i = 0; i < 8; ++i, coord += velocity) {
		if (coord.s > 1.0 || coord.t > 1.0 || coord.s < 0.0 || coord.t < 0.0) {
			break;
		}
		if (texture2D(gaux1, coord).r < 0.5) {
			color += texture2D(composite, coord);
			++samples;
		}
	}
*/
		
	
#endif



/*
#ifdef BLOOM
	color = color * 0.8;
	color += addBloom(color, texcoord.st);
#endif
*/

#ifdef SSAO
  float lum = dot(color.rgb, vec3(1.0));
  vec3 luminance = vec3(lum);
  //vec3 ssaodark = getSSAOFactor();
  color.rgb *= mix(getSSAOFactor(), vec3(1.0), luminance * SSAO_LUMINANCE);
  //color.r = color.r * (pow(ssaodark.r, 1.5));
  //color.g = color.g * (pow(ssaodark.g, 1.5));
  //color.b = color.b * (pow(ssaodark.b, 1.5));
#endif

#ifdef BLOOM

	color = color * 0.8;
	
	float radius = 0.002;
	float blm_amount = 0.02*BLOOM_AMOUNT;
	float sc = 20.0;
	
	int i = 0;
	int samples = 1;
	
	vec4 clr = vec4(0.0);
	
	for (i = -10; i < 10; i++) {
	clr += texture2D(composite, texcoord.st + (vec2(i,i))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(i,-i))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(-i,i))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(-i,-i))*radius)*sc;
	
	clr += texture2D(composite, texcoord.st + (vec2(0.0,i))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-i))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(-i,0.0))*radius)*sc;
	clr += texture2D(composite, texcoord.st + (vec2(i,0.0))*radius)*sc;
	
	++samples;
	sc = sc - 1.0;
	}
	
	clr = (clr/8.0)/samples;
	
	color += clr*blm_amount;

#endif

#ifdef GLARE

	color = color * 0.8;
	
	float radius = 0.002*GLARE_RANGE;
	float radiusv = 0.002;
	float bloomintensity = 0.1*GLARE_AMOUNT;

	vec4 clr = vec4(0.0);
	
	clr += texture2D(composite, texcoord.st);
	
	//horizontal (70 taps)
	
	clr += texture2D(composite, texcoord.st + (vec2(19.0,0.0))*radius)*1.0;
	clr += texture2D(composite, texcoord.st + (vec2(18.0,0.0))*radius)*2.0;
	clr += texture2D(composite, texcoord.st + (vec2(17.0,0.0))*radius)*3.0;
	clr += texture2D(composite, texcoord.st + (vec2(16.0,0.0))*radius)*4.0;
	clr += texture2D(composite, texcoord.st + (vec2(15.0,0.0))*radius)*5.0;
	clr += texture2D(composite, texcoord.st + (vec2(14.0,0.0))*radius)*6.0;
	clr += texture2D(composite, texcoord.st + (vec2(13.0,0.0))*radius)*7.0;
	clr += texture2D(composite, texcoord.st + (vec2(12.0,0.0))*radius)*8.0;
	clr += texture2D(composite, texcoord.st + (vec2(11.0,0.0))*radius)*9.0;
	clr += texture2D(composite, texcoord.st + (vec2(10.0,0.0))*radius)*10.0;
	clr += texture2D(composite, texcoord.st + (vec2(9.0,0.0))*radius)*11.0;
	clr += texture2D(composite, texcoord.st + (vec2(8.0,0.0))*radius)*12.0;
	clr += texture2D(composite, texcoord.st + (vec2(7.0,0.0))*radius)*13.0;
	clr += texture2D(composite, texcoord.st + (vec2(6.0,0.0))*radius)*14.0;
	clr += texture2D(composite, texcoord.st + (vec2(5.0,0.0))*radius)*15.0;
	clr += texture2D(composite, texcoord.st + (vec2(4.0,0.0))*radius)*16.0;
	clr += texture2D(composite, texcoord.st + (vec2(3.0,0.0))*radius)*17.0;
	clr += texture2D(composite, texcoord.st + (vec2(2.0,0.0))*radius)*18.0;
	clr += texture2D(composite, texcoord.st + (vec2(1.0,0.0))*radius)*19.0;
	
		clr += texture2D(composite, texcoord.st + (vec2(0.0,0.0))*radius)*20.0;
		
	clr += texture2D(composite, texcoord.st + (vec2(-1.0,0.0))*radius)*19.0;
	clr += texture2D(composite, texcoord.st + (vec2(-2.0,0.0))*radius)*18.0;
	clr += texture2D(composite, texcoord.st + (vec2(-3.0,0.0))*radius)*17.0;
	clr += texture2D(composite, texcoord.st + (vec2(-4.0,0.0))*radius)*16.0;
	clr += texture2D(composite, texcoord.st + (vec2(-5.0,0.0))*radius)*15.0;
	clr += texture2D(composite, texcoord.st + (vec2(-6.0,0.0))*radius)*14.0;
	clr += texture2D(composite, texcoord.st + (vec2(-7.0,0.0))*radius)*13.0;
	clr += texture2D(composite, texcoord.st + (vec2(-8.0,0.0))*radius)*12.0;
	clr += texture2D(composite, texcoord.st + (vec2(-9.0,0.0))*radius)*11.0;
	clr += texture2D(composite, texcoord.st + (vec2(-10.0,0.0))*radius)*10.0;
	clr += texture2D(composite, texcoord.st + (vec2(-11.0,0.0))*radius)*9.0;
	clr += texture2D(composite, texcoord.st + (vec2(-12.0,0.0))*radius)*8.0;
	clr += texture2D(composite, texcoord.st + (vec2(-13.0,0.0))*radius)*7.0;
	clr += texture2D(composite, texcoord.st + (vec2(-14.0,0.0))*radius)*6.0;
	clr += texture2D(composite, texcoord.st + (vec2(-15.0,0.0))*radius)*5.0;
	clr += texture2D(composite, texcoord.st + (vec2(-16.0,0.0))*radius)*4.0;
	clr += texture2D(composite, texcoord.st + (vec2(-17.0,0.0))*radius)*3.0;
	clr += texture2D(composite, texcoord.st + (vec2(-18.0,0.0))*radius)*2.0;
	clr += texture2D(composite, texcoord.st + (vec2(-19.0,0.0))*radius)*1.0;
	
	//vertical
	clr += texture2D(composite, texcoord.st + (vec2(0.0,19.0))*radius)*1.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,18.0))*radius)*2.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,17.0))*radius)*3.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,16.0))*radius)*4.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,15.0))*radius)*5.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,14.0))*radius)*6.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,13.0))*radius)*7.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,12.0))*radius)*8.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,11.0))*radius)*9.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,10.0))*radius)*10.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,9.0))*radius)*11.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,8.0))*radius)*12.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,7.0))*radius)*13.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,6.0))*radius)*14.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,5.0))*radius)*15.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,4.0))*radius)*16.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,3.0))*radius)*17.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,2.0))*radius)*18.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,1.0))*radius)*19.0;
	
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-19.0))*radius)*1.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-18.0))*radius)*2.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-17.0))*radius)*3.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-16.0))*radius)*4.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-15.0))*radius)*5.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-14.0))*radius)*6.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-13.0))*radius)*7.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-12.0))*radius)*8.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-11.0))*radius)*9.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-10.0))*radius)*10.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-9.0))*radius)*11.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-8.0))*radius)*12.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-7.0))*radius)*13.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-6.0))*radius)*14.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-5.0))*radius)*15.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-4.0))*radius)*16.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-3.0))*radius)*17.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-2.0))*radius)*18.0;
	clr += texture2D(composite, texcoord.st + (vec2(0.0,-1.0))*radius)*19.0;

	clr = (clr/77.0)/5.0;
	clr.r = pow(clr.r, 1.2)*1.6 - (clr.g + clr.b)*0.6;
	clr.g = pow(clr.g, 1.2)*1.6 - (clr.r + clr.b)*0.6;
	clr.b = pow(clr.b, 1.2)*1.6 - (clr.r + clr.g)*0.6;
	
	clr = clamp((clr), 0.0, 1.0);
	
	color.r = color.r + (clr.r*1.5)*bloomintensity;
	color.g = color.g + (clr.g*1.5)*bloomintensity;
	color.b = color.b + (clr.b*4.0)*bloomintensity;
	color = max(color, 0.0);
	//color = color*1.05 - 0.05;
	
	

#endif


#ifdef CEL_SHADING
	color.rgb *= (getCellShaderFactor(texcoord.st));
#endif

#ifdef HDR

float avgclr = realcolor(vec2(0.5, 0.5));
      avgclr += realcolor(vec2(0.2, 0.2));
	  avgclr += realcolor(vec2(0.2, -0.2));
	  avgclr += realcolor(vec2(-0.2, 0.2));
	  avgclr += realcolor(vec2(-0.2, -0.2));
	  avgclr = clamp(avgclr/5, 0.0, 0.8)*1.2;
	  
color.r = color.r + color.r*((1 - avgclr)*2);
color.g = color.g + color.g*((1 - avgclr)*2);
color.b = color.b + color.b*((1 - avgclr)*2);
	  


#endif;

#ifdef CROSSPROCESS
	//pre-gain
	color.r = color.r * (BRIGHTMULT + 0.2);
	color.g = color.g * (BRIGHTMULT + 0.2);
	color.b = color.b * (BRIGHTMULT + 0.2);
	
	//compensate for low-light artifacts
	color = color+0.025;

	//calculate double curve
	float dbr = -color.r + 1.4;
	float dbg = -color.g + 1.4;
	float dbb = -color.b + 1.4;
	
	//fade between simple gamma up curve and double curve
	float pr = mix(dbr, 0.55, 0.7);
	float pg = mix(dbg, 0.55, 0.7);
	float pb = mix(dbb, 0.85, 0.7);
	
	color.r = pow((color.r * 0.99 - 0.02), pr);
	color.g = pow((color.g * 0.99 - 0.015), pg);
	color.b = pow((color.b * 0.7 + 0.04), pb);
#endif

#ifdef HIGHDESATURATE

	//desaturate technique (choose one)

	//average
	float rgb = max(color.r, max(color.g, color.b))/2 + min(color.r, min(color.g, color.b))/2;

	//adjust black and white image to be brighter
	float bw = pow(rgb, 0.7);

	//mix between per-channel analysis and average analysis
	float rgbr = mix(rgb, color.r, 0.7);
	float rgbg = mix(rgb, color.g, 0.7);
	float rgbb = mix(rgb, color.b, 0.7);

	//calculate crossfade based on lum
	float mixfactorr = max(0.0, (rgbr*3 - 2));
	float mixfactorg = max(0.0, (rgbg*3 - 2));
	float mixfactorb = max(0.0, (rgbb*3 - 2));

	//crossfade between saturated and desaturated image
	float mixr = mix(color.r, bw, mixfactorr);
	float mixg = mix(color.g, bw, mixfactorg);
	float mixb = mix(color.b, bw, mixfactorb);

	//adjust level of desaturation
	color.r = clamp((mix(mixr, color.r, 0.2)), 0.0, 1.0);
	color.g = clamp((mix(mixg, color.g, 0.2)), 0.0, 1.0);
	color.b = clamp((mix(mixb, color.b, 0.2)), 0.0, 1.0);
	

	//hold color values for color boost
	//vec4 hld = color;

	
	//Color boosting
	color.r = (color.r)*(COLOR_BOOST + 1.0) + (color.g + color.b)*(-COLOR_BOOST);
	color.g = (color.g)*(COLOR_BOOST + 1.0) + (color.r + color.b)*(-COLOR_BOOST);
	color.b = (color.b)*(COLOR_BOOST + 1.0) + (color.r + color.g)*(-COLOR_BOOST);
	
	//color.r = mix(((color.r)*(COLOR_BOOST + 1.0) + (hld.g + hld.b)*(-COLOR_BOOST)), hld.r, (max(((1-rgb)*2 - 1), 0.0)));
	//color.g = mix(((color.g)*(COLOR_BOOST + 1.0) + (hld.r + hld.b)*(-COLOR_BOOST)), hld.g, (max(((1-rgb)*2 - 1), 0.0)));
	//color.b = mix(((color.b)*(COLOR_BOOST + 1.0) + (hld.r + hld.g)*(-COLOR_BOOST)), hld.b, (max(((1-rgb)*2 - 1), 0.0)));
	
	//undo artifact compensation
	color = max(((color*1.13)-0.03), 0.0);
	color = color*1.02 - 0.02;

	
	
#endif
	gl_FragColor = color;
	
// End of Main. -----------------
}
