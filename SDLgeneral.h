
#define mapSize 18
#define N 18

//SDL initial and close
SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;
SDL_Texture* walltexture = NULL;
SDL_Texture* enemytexture = NULL;
SDL_Texture* presstexture = NULL;
SDL_Texture* cursortexture = NULL;
SDL_Texture* skeletontexture = NULL;
SDL_Texture* skeleton_redtexture = NULL;
SDL_Texture* porttexture = NULL;
SDL_Texture* ponpon_1_texture = NULL;
SDL_Texture* ponpon_2_texture = NULL;
SDL_Texture* waytexture = NULL;


//*** MUSIC ***
/*
Mix_Music* entMusic = NULL;
Mix_Music* menuMusic = NULL;
Mix_Music* jesusMusic = NULL;
Mix_Music* failMusic = NULL;
Mix_Music* survivalMusic = NULL;
Mix_Music* pauseMusic = NULL;
Mix_Music* successMusic = NULL;
Mix_Music* alertMusic = NULL;
*/
//*** CHUNK ***
Mix_Chunk* entMusic = NULL; //channel 1
Mix_Chunk* menuMusic = NULL; //channel 2
Mix_Chunk* jesusMusic = NULL; //channel 3
Mix_Chunk* failMusic = NULL; //channel 4
Mix_Chunk* survivalMusic = NULL; //channel 5
Mix_Chunk* pauseMusic = NULL; //channel 6
Mix_Chunk* successMusic = NULL; //channel 7
Mix_Chunk* alertMusic = NULL; //channel 8

Mix_Chunk* transMusic = NULL; //channel 9
Mix_Chunk* clickMusic = NULL; //channel 10
Mix_Chunk* timerMusic = NULL; //channel 11

bool survivalbgm;
bool timerbgm;
bool jesusbgm;
bool failbgm;
bool pausebgm;
bool successbgm;
bool alertbgm;
bool transbgm;
//click has no need

//screen dimension
const int WIDTH = 1200;
const int HEIGHT = 750;
const int sWIDTH = 200;
const int sHEIGHT = 250;

const int FirstW = WIDTH - 200;
const int FirstH = HEIGHT;

int initSDL()
{
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0)
	{
		printf("SDL_Init failed: %s\n", SDL_GetError());
		return 1;
	}
	// Create window
	// SDL_WINDOWPOS_UNDEFINED: Used to indicate that you don't care what the window position is.
	window = SDL_CreateWindow("BATTLE ROYALE 3D", 50, 50, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	if (window == NULL)
	{
		printf("SDL_CreateWindow failed: %s\n", SDL_GetError());
		SDL_Quit();
		return 2;
	}
	int imgFlags = IMG_INIT_PNG;
	if (!(IMG_Init(imgFlags) & imgFlags))
	{
		printf("SDL_image failed: %s\n", IMG_GetError());
		return 3;
	}

	if (TTF_Init() == -1)
	{
		// https://www.libsdl.org/projects/SDL_ttf/docs/SDL_ttf_frame.html
		printf("SDL_ttf could not initialize! SDL_ttf Error: %s\n", TTF_GetError());
		return 4;
	}
	// Initialize SDL_mixer
	if (Mix_OpenAudio(22050, MIX_DEFAULT_FORMAT, 8, 2048) < 0)
	{
		printf("SDL_mixer could not initialize! SDL_mixer Error: %s\n", Mix_GetError());
		return 5;
	}
	Mix_AllocateChannels(12);
	// Create renderer
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	if (renderer == NULL)
	{
		SDL_DestroyWindow(window);
		printf("SDL_CreateRenderer failed: %s\n", SDL_GetError());
		SDL_Quit();
		return 6;
	}
	return 0;
}
void closeSDL()
{
	/*
	//FREE MUSIC
	Mix_FreeMusic(entMusic);
	Mix_FreeMusic(menuMusic);
	Mix_FreeMusic(jesusMusic);
	Mix_FreeMusic(failMusic);
	Mix_FreeMusic(survivalMusic);
	Mix_FreeMusic(pauseMusic);
	Mix_FreeMusic(successMusic);
	Mix_FreeMusic(alertMusic);
	*/

	//FREE CHUNK
	Mix_FreeChunk(entMusic);
	Mix_FreeChunk(menuMusic);
	Mix_FreeChunk(jesusMusic);
	Mix_FreeChunk(failMusic);
	Mix_FreeChunk(survivalMusic);
	Mix_FreeChunk(pauseMusic);
	Mix_FreeChunk(successMusic);
	Mix_FreeChunk(alertMusic);

	Mix_FreeChunk(transMusic);
	Mix_FreeChunk(clickMusic);
	Mix_FreeChunk(timerMusic);

	// Destroy renderer
	// Destroy window
	// Quit Image subsystem
	// Quit SDL subsystems
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	// Shutdown subsystem
	Mix_Quit();
	TTF_Quit();
	SDL_Quit();
}

bool loadAudio()
{
	/*
	// Load music
	entMusic = Mix_LoadMUS("./bgm/Home.mp3");
	menuMusic = Mix_LoadMUS("./bgm/Menu.mp3");
	jesusMusic = Mix_LoadMUS("./bgm/God.mp3");
	survivalMusic = Mix_LoadMUS("./bgm/Man_Survival.mp3");
	pauseMusic = Mix_LoadMUS("./bgm/Pause.mp3");
	successMusic = Mix_LoadMUS("./bgm/Success.mp3");
	failMusic = Mix_LoadMUS("./bgm/Fail.mp3");
	alertMusic = Mix_LoadMUS("./bgm/Alert.mp3");
	if (entMusic == NULL || menuMusic == NULL || alertMusic == NULL || jesusMusic == NULL || survivalMusic == NULL || pauseMusic == NULL || successMusic == NULL || failMusic == NULL)
	{
		printf("1. Failed to load music! SDL_mixer Error: %s\n", Mix_GetError());
		return false;
	}
	*/

	// Load sound effects
	entMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Home.wav");
	menuMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Menu.wav");
	jesusMusic = Mix_LoadWAV("./bgm/Bgm_War Series/God.wav");
	survivalMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Man_Survival.wav");
	pauseMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Pause.wav");
	successMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Success.wav");
	failMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Fail.wav");
	alertMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Alert.wav");

	clickMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Click.wav");
	timerMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Man_Timer.wav");
	transMusic = Mix_LoadWAV("./bgm/Bgm_War Series/Loading.wav");
	if (entMusic == NULL || menuMusic == NULL || alertMusic == NULL || jesusMusic == NULL || survivalMusic == NULL || pauseMusic == NULL || successMusic == NULL || failMusic == NULL || clickMusic == NULL || timerMusic == NULL || transMusic == NULL)
	{
		printf("Failed to load sound effect! SDL_mixer Error: %s\n", Mix_GetError());
		return false;
	}

}
//texture
const int SOLID = 100;
const int SHADED = 101;
const int BLENDED = 102;

struct TextData
{
	SDL_Texture* texture;
	int width;
	int height;
};
struct ImageData
{
	SDL_Texture* texture;
	int width;
	int height;
	int num;
	int wn;
	int hn;
};
struct ImgPos {
	SDL_Rect src;
	SDL_Rect dst;
};
struct text {
	int x, y;
	int w, h;
	int r, g, b, a;
	bool hover = false;
	bool click = false;
};
enum MouseState
{
	NONE = 0,
	OUT = 1, // Mouse out
	IN_LB_SC = 2,  // Inside, Left Button, Single Click
	IN_RB_SC = 3,  // Inside, RIGHT Button, Single Click
	IN_LB_DC = 4,  // Inside, Left Button, Double Click
	IN_RB_DC = 5,  // Inside, RIGHT Button, Double Click
	IN_LB_PR = 6,  // Inside, Left Button, Press
	IN_RB_PR = 7,  // Inside, Left Button, Press
	IN_WU = 8,  // Inside, Wheel UP
	IN_WD = 9,  // Inside, Wheel DOWN
	HOVER = 10, // Mouse hover
	IN_LB_PR_HOVER = 11
};


TextData loadTextTexture(const char* str, const char* fontPath, int fontSize, Uint8 fr, Uint8 fg, Uint8 fb, int textType, Uint8 br, Uint8 bg, Uint8 bb)
{
	TTF_Font* ttfFont = NULL;
	ttfFont = TTF_OpenFont(fontPath, fontSize);
	if (ttfFont == NULL)
	{
		printf("Failed to load lazy font! SDL_ttf Error: %s\n", TTF_GetError());
	}
	SDL_Color textFgColor = { fr, fg, fb }, textBgColor = { br, bg, bb };
	SDL_Surface* textSurface = NULL;
	if (textType == SOLID)
	{
		textSurface = TTF_RenderText_Solid(ttfFont, str, textFgColor);
	}
	else if (textType == SHADED)
	{
		textSurface = TTF_RenderText_Shaded(ttfFont, str, textFgColor, textBgColor);
	}
	else if (textType == BLENDED)
	{
		textSurface = TTF_RenderText_Blended(ttfFont, str, textFgColor);
	}

	TTF_CloseFont(ttfFont);

	if (textSurface == NULL)
	{
		printf("Unable to render text surface! SDL_ttf Error: %s\n", TTF_GetError());
		return { NULL };
	}
	else
	{
		TextData text;
		text.texture = SDL_CreateTextureFromSurface(renderer, textSurface);
		if (text.texture == NULL)
		{
			printf("SDL_CreateTextureFromSurface failed: %s\n", SDL_GetError());
		}
		text.width = textSurface->w;
		text.height = textSurface->h;
		SDL_FreeSurface(textSurface);

		return text;
	}
}
ImageData loadImgTexture(char* path, int num, int hn, int wn, bool ckEnable, Uint8 r, Uint8 g, Uint8 b)
{
	//Load image at specified path
	SDL_Surface* loadedSurface = IMG_Load(path);

	if (loadedSurface == NULL)
	{
		printf("IMG_Load failed: %s\n", IMG_GetError());
		return { NULL };
	}
	else
	{
		ImageData img;

		// Set the color key (transparent pixel) in a surface.		
		SDL_SetColorKey(loadedSurface, ckEnable, SDL_MapRGB(loadedSurface->format, r, g, b));

		// Create texture from surface pixels
		img.texture = SDL_CreateTextureFromSurface(renderer, loadedSurface);
		if (img.texture == NULL)
		{
			printf("SDL_CreateTextureFromSurface failed: %s\n", SDL_GetError());
		}

		//Get image dimensions and information
		img.width = loadedSurface->w;
		img.height = loadedSurface->h;
		img.num = num;
		img.wn = wn;
		img.hn = hn;

		// Get rid of old loaded surface
		SDL_FreeSurface(loadedSurface);

		//return newTexture;
		return img;
	}
}
void imgRender(SDL_Renderer* renderer, ImageData img, int posX, int posY)
{
	SDL_Rect r;
	r.x = posX;
	r.y = posY;
	r.w = img.width;
	r.h = img.height;
	SDL_RenderCopy(renderer, img.texture, NULL, &r);
}
int textRender(SDL_Renderer* renderer, TextData text, int posX, int posY, int cx, int cy, double angle, SDL_RendererFlip flip, int alpha)
{
	SDL_Rect r;
	r.x = posX;
	r.y = posY;
	r.w = text.width;
	r.h = text.height;

	if (SDL_SetTextureBlendMode(text.texture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}

	if (SDL_SetTextureAlphaMod(text.texture, alpha) == -1)
	{
		printf("SDL_SetTextureAlphaMod failed: %s\n", SDL_GetError());
		return -1;
	}

	SDL_Point center = { cx, cy };

	SDL_RenderCopyEx(renderer, text.texture, NULL, &r, angle, &center, flip);
	return 1;
}

//****************** ENTRANCE *******************
bool entFlag = true; //back to entrance
bool entbgm = true;
Uint32 TitleDrop(Uint32 interval, void* param) {
	int* par = (int*)param;
	par[1] += 3;
	if (par[1] > HEIGHT / 2 - par[2] / 2 - 40) {

		return 0;
	}
	return interval;
}


//************** MENU *********************
bool menuFlag = false;
bool menubgm = false;
text t_HOME, t_menuTitle, t_MODE, t_WORLD, t_MUSIC, t_QUIT, c_INVERSION, c_TIME, c_survival, c_simple, c_classic, c_on, c_off, t_PLAY;

//************** TRANSFORM *********************
bool transFlag = false;
int duration = 10000; //first: 10 sec, after: 3 sec
clock_t trans_start = 0;

//************** JESUS *********************
bool backgroundDisplay = true;
bool jesusDisplay = true;
bool jesusFlag = false; //enter jesus perspective
bool recentFace = false;

//background
int starX[100] = { 0 };
int starY[100] = { 0 };

char way[100] = "images/way.png";
ImageData way_img;

char ponpon_Path1[100] = "images/ponpon_1.png";
char ponpon_Path2[100] = "images/ponpon_2.png";
ImageData ponpon_1;
ImageData ponpon_2;

//****************** FIRST *******************
float depthBuffer[sWIDTH];

float shadowC[sHEIGHT][sWIDTH] = { 0 };
float shadowW[sHEIGHT][sWIDTH] = { 0 };
float shadowG[sHEIGHT][sWIDTH] = { 0 };
float shadowE[sHEIGHT][sWIDTH] = { 0 };
float shadowP[sHEIGHT][sWIDTH] = { 0 };

//counting time when game start
int startTime;
int now;
int pauseTime;

bool firstDisplay = true;
bool firstFlag = false; //enter first-perspective
bool end = false;
bool H_pause = false;
bool pauseSurface = false;//pause surface
bool approach = false;
bool quickFlag = false;
bool reveal = true;
bool alertflag = true;

text t_f_time, t_f_pos, t_f_title, t_f_press;

//char wallimg[100] = "images/woodwall.png";
//char wallimg[100] = "images/corkwall.png";
//char wallimg[100] = "images/stonewall.png";
//char wallimg[100] = "images/square.png";
char wallimg[100] = "images/brickwall.png";
ImageData wallpaper;
char enemyimg[4][100] = { "images/Enemy1_1.png", "images/Enemy2_1.png", "images/Enemy3.png", "images/Enemy4.png" };
//char enemyimg[100] = "images/Enemy2.png";
//char enemyimg[100] = "images/smallghost.png";
ImageData enemypaper;
char skeletonimg[100] = "images/skeleton.png";
ImageData skeletonpaper;
char skeleton_redimg[100] = "images/skeleton_red.png";
ImageData skeleton_redpaper;
char portimg[100] = "images/light3.png";
ImageData portpaper;


//************** PAUSE SURFACE *********************
text t_continue, t_back;

//************** FIRST && END *********************
bool succeed = false, fail = false, timeEND = false;
text t_hash, t_AND, t_zero;
text c_BackToMenu, c_TryAgain;
const int game_over_SizeW = 64;
const int game_over_SizeH = 22;
int timeLimit = 120;

char game_over[game_over_SizeH][game_over_SizeW] =
{
	{'&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&'},
	{'&', '&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&', '&'},
	{'&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&'},
	{'&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '#', '#', '#', '#', '#', '0', '#', '#', '0', '0', '#', '0', '0', '0', '#', '0', '0', '0', '0', '#', '0', '0', '0', '#', '#', '#', '0', '0', '0', '#', '0', '0', '0', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '#', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '#', '0', '0', '0', '0', '#', '0', '#', '0', '0', '0', '#', '0', '0', '#', '#', '0', '0', '#', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '#', '#', '#', '0', '0', '#', '#', '#', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '#', '0', '#', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '#', '0', '#', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '#', '0', '0', '0', '0', '#', '0', '0', '0', '0', '#', '#', '#', '#', '#', '0', '#', '0', '#', '0', '#', '#', '#', '#', '#', '0', '0', '#', '0', '0', '#', '0', '0', '#', '#', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '#', '0', '0', '#', '0', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&'},
	{'&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&'},
	{'&', '&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&', '&'},
	{'&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&'},
	{'&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&'},
	{'&', '&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&', '&'},
	{'&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&'},
	{'&', '0', '0', '0', '0', '#', '#', '#', '0', '0', '0', '0', '#', '0', '0', '0', '#', '#', '#', '0', '#', '0', '0', '#', '0', '0', '0', '#', '#', '#', '#', '#', '#', '#', '#', '#', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '#', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '#', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '0', '0', '#', '#', '0', '#', '#', '0', '#', '0', '0', '0', '#', '#', '0', '0', '#', '0', '#', '0', '0', '#', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '#', '#', '#', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '0', '#', '#', '0', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '#', '0', '#', '#', '#', '0', '#', '0', '#', '0', '#', '0', '#', '0', '0', '#', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '#', '0', '0', '#', '0', '#', '#', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '0', '0', '0', '0', '#', '0', '0', '#', '0', '0', '#', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '0', '0', '0', '#', '0', '0', '#', '#', '0', '#', '0', '0', '#', '0', '0', '0', '0', '&'},
	{'&', '0', '0', '0', '0', '#', '#', '#', '0', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '0', '#', '0', '0', '#', '#', '0', '0', '0', '0', '#', '0', '0', '#', '#', '#', '#', '0', '0', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '0', '#', '0', '0', '0', '#', '0', '#', '#', '#', '#', '0', '0', '0', '0', '&'},
	{'&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&'},
	{'&', '&', '&', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '&', '&', '&'},
	{'&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&', '&'}
};

/*
flag record:
	1. entFlag + entbgm => entrance interface
	2. menuFlag + menubgm => menu interface
	3. jesusFlag + jesusDisplay => jesus interface
	4. firstFlag + firstDisplay => first interface
	5. succeed
	6. fail
*/

Uint32 approachflash(Uint32 interval, void* param) {


	return interval;
}

Uint32 background(Uint32 interval, void* param) {

	srand(time(NULL));
	for (int i = 0; i < 100; i++) {
		starX[i] = rand() % WIDTH;
		starY[i] = rand() % HEIGHT;
	}

	return interval;
}


//Handle Event
double fFov = 90.0;

//****************** FIRST *******************
clock_t t;
float elapsedTime;

float userX = 5.0f;		                // Player Start Position
float userY = 5.0f;
float fov = M_PI / 3.f;             	// Field of View
float userR = M_PI;			            // Player Start Rotation
float depth = (float)mapSize;		    // Maximum rendering distance cuz mapSize is 10*10
float speed = 30.0f;			        // Walking Speed
float shiftR = 0, shiftX = 0, shiftY = 0;

float rayAngle, rayStepSize = 0.01f, dToWall, dToEnemy, dToPort;
bool hitWall, boundary, hitEnemy, hitPort;
float eyeX, eyeY;
int testX, testY;
float WsampleX = 0.0f;
float EsampleX = 0.0f;
float PsampleX = 0.0f;
ImgPos wallPos[sHEIGHT][sWIDTH], EImgPos[sHEIGHT][sWIDTH], PImgPos[sHEIGHT][sWIDTH];

float Espeed = 8.0f;
float enemyX = 11.0f;		// enemy Start Position
float enemyY = 9.0f;
float enemyR = 0.0f;		// enemy Start Rotation
float Efov = M_PI / 4.0f;	// Field of View
float EshiftR = 0, EshiftX = 0, EshiftY = 0;
float ErayAngle, ErayStepSize = 0.01f, EdToWall;
bool EhitWall, Eboundary;
float EeyeX, EeyeY, dToUser;
int EtestX, EtestY;
bool hitUser;

int portX, portY;

int mouse_lastX = WIDTH / 2, mouse_lastY = HEIGHT / 2;



void mouseHandleEvent(SDL_Event* e, MouseState* mouseState, int* x, int* y)
{

	if (e->type == SDL_MOUSEMOTION || e->type == SDL_MOUSEBUTTONDOWN || e->type == SDL_MOUSEBUTTONUP || e->type == SDL_MOUSEWHEEL)
	{
		SDL_GetMouseState(x, y);

		static int lasttime = SDL_GetTicks();
		static int curtime = SDL_GetTicks();
		int timediv;

		lasttime = curtime;
		curtime = SDL_GetTicks();
		timediv = curtime - lasttime;

		switch (e->type)
		{
		case SDL_MOUSEBUTTONDOWN:
			break;

		case SDL_MOUSEBUTTONUP:
			//double clicks
			if (e->button.button == SDL_BUTTON_LEFT && e->button.clicks == 2 && timediv < 150)
			{
				*mouseState = IN_LB_DC;
			}
			else if (e->button.button == SDL_BUTTON_RIGHT && e->button.clicks == 2 && timediv < 150)
			{
				*mouseState = IN_RB_DC;
			}

			//single click
			else if (e->button.button == SDL_BUTTON_LEFT && e->button.clicks == 1 && timediv < 800 && timediv > 5)
			{
				//menu
				if (menuFlag) {
					Mix_PlayChannel(10, clickMusic, 0);
					if (*x > t_QUIT.x && *x<(t_QUIT.x + t_QUIT.w) && *y>t_QUIT.y && *y < (t_QUIT.y + t_QUIT.h)) {
						exit(1);
					}
					if (*x > t_HOME.x && *x<(t_HOME.x + t_HOME.w) && *y>t_HOME.y && *y < (t_HOME.y + t_HOME.h)) {
						entFlag = true;
						menuFlag = false;
						Mix_HaltChannel(2);
						Mix_PlayChannel(1, entMusic, -1);
					}
					if (*x > c_TIME.x && *x<(c_TIME.x + c_TIME.w) && *y>c_TIME.y && *y < (c_TIME.y + c_TIME.h)) {
						c_TIME.click = true;
						c_survival.click = false;
						c_INVERSION.click = false;
					}
					if (*x > c_survival.x && *x<(c_survival.x + c_survival.w) && *y>c_survival.y && *y < (c_survival.y + c_survival.h)) {
						c_survival.click = true;
						c_TIME.click = false;
						c_INVERSION.click = false;
					}
					if (*x > c_simple.x && *x<(c_simple.x + c_simple.w) && *y>c_simple.y && *y < (c_simple.y + c_simple.h)) {
						c_simple.click = true;
						c_classic.click = false;
					}
					if (*x > c_classic.x && *x<(c_classic.x + c_classic.w) && *y>c_classic.y && *y < (c_classic.y + c_classic.h)) {
						c_classic.click = true;
						c_simple.click = false;
					}
					if (*x > c_on.x && *x<(c_on.x + c_on.w) && *y>c_on.y && *y < (c_on.y + c_on.h)) {
						c_on.click = true;
						c_off.click = false;
					}
					if (*x > c_off.x && *x<(c_off.x + c_off.w) && *y>c_off.y && *y < (c_off.y + c_off.h)) {
						c_off.click = true;
						c_on.click = false;
					}
					if (*x > c_INVERSION.x && *x<(c_INVERSION.x + c_INVERSION.w) && *y>c_INVERSION.y && *y < (c_INVERSION.y + c_INVERSION.h)) {
						c_INVERSION.click = true;
						c_survival.click = false;
						c_TIME.click = false;
					}
					//PRESS PLAY
					if (*x > t_PLAY.x && *x<(t_PLAY.x + t_PLAY.w) && *y>t_PLAY.y && *y < (t_PLAY.y + t_PLAY.h)) {
						transFlag = true;
						transbgm = true;
						entFlag = false;
						menuFlag = false;
						Mix_HaltChannel(2);
						if (transbgm) {
							Mix_PlayChannel(9, transMusic, -1);
							transbgm = false;
						}
					}
				}

				//first sight
				if (firstFlag && !(fail || succeed)) {
					if (sqrt(pow(*x - 100, 2) + powf((*y - 620), 2)) <= 50) {
						Mix_PlayChannel(10, clickMusic, 0);
						if (!pauseSurface) {
							//if (c_survival.click) Mix_HaltChannel(5);
							//else Mix_HaltChannel(11);
							Mix_FadeOutChannel(-1, 1500);
							Mix_PlayChannel(6, pauseMusic, -1);
							pauseSurface = true;
							pauseTime = clock();
						}
						else {
							Mix_HaltChannel(6);
							if (c_survival.click) {
								//Mix_Resume(5);
								//if (approach)Mix_Resume(8);
								Mix_FadeInChannelTimed(5, survivalMusic, -1, 1500, -1);
								if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
								//Mix_PlayChannel(5, survivalMusic, -1);
								//if (approach)Mix_PlayChannel(8, alertMusic, -1);
							}
							else {
								//Mix_Resume(11);
								//if (approach)Mix_Resume(8);
								Mix_FadeInChannelTimed(11, timerMusic, -1, 1500, -1);
								if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
								//Mix_PlayChannel(11, timerMusic, -1);
								//if (approach)Mix_PlayChannel(8, alertMusic, -1);
							}
							pauseSurface = false;
							startTime += (clock() - pauseTime);
						}
					}
				}

				//jesus
				else if ((jesusFlag) && !(fail || succeed)) {
					if (sqrt(pow(*x - 100, 2) + powf((*y - 620), 2)) <= 50) {
						Mix_PlayChannel(10, clickMusic, 0);
						if (!pauseSurface) {
							Mix_FadeOutChannel(3, 1500);
							//Mix_HaltChannel(3);
							Mix_PlayChannel(6, pauseMusic, -1);
							pauseSurface = true;
						}
						else {
							Mix_HaltChannel(6);
							Mix_FadeInChannelTimed(3, jesusMusic, -1, 1500, -1);
							if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
							//Mix_PlayChannel(3, jesusMusic, -1);
							pauseSurface = false;
						}
					}
				}

				if ((fail || succeed)) {
					//play again
					if (*x > (t_AND.x + 12) && *x<(t_AND.x + (game_over_SizeW - 1) * 12) && *y>(t_AND.y + 12) && *y < (t_AND.y + 12 * (game_over_SizeH / 2 - 1))) {
						Mix_PlayChannel(10, clickMusic, 0);
						firstFlag = true;
						firstDisplay = true;
						jesusDisplay = false;
						fail = false;
						succeed = false;
						/*if (fail) Mix_HaltChannel(4);
						else Mix_HaltChannel(7);*/
						Mix_HaltChannel(-1);
						;						if (c_survival.click) Mix_PlayChannel(5, survivalMusic, -1);
						else Mix_PlayChannel(11, timerMusic, -1);
					}
					// back to menu
					else if (*x > (t_AND.x + 12) && *x<(t_AND.x + (game_over_SizeW - 1) * 12) && *y>(t_AND.y + 12 * (game_over_SizeH / 2 + 1)) && *y < (t_AND.y + 12 * (game_over_SizeH - 1))) {
						Mix_PlayChannel(10, clickMusic, 0);
						menuFlag = true;
						firstDisplay = false;
						firstFlag = false;
						/*if (fail) Mix_HaltChannel(4);
						else Mix_HaltChannel(7);*/
						Mix_HaltChannel(-1);
						Mix_PlayChannel(2, menuMusic, -1);
						fail = false;
						succeed = false;
						c_TIME.click = true;
						c_survival.click = false;
						c_INVERSION.click = false;
						c_simple.click = true;
						c_classic.click = false;
					}
				}

				//pause surface
				if (pauseSurface) {
					if (*x > t_continue.x && *x<(t_continue.x + t_continue.w) && *y>t_continue.y && *y < (t_continue.y + t_continue.h)) {
						Mix_PlayChannel(10, clickMusic, 0);
						Mix_HaltChannel(-1);
						if (jesusFlag) {
							Mix_FadeInChannelTimed(3, jesusMusic, -1, 1500, -1);
							if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
							//Mix_PlayChannel(3, jesusMusic, -1);
						}
						else {
							if (c_survival.click) {
								Mix_FadeInChannelTimed(5, survivalMusic, -1, 1500, -1);
								if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
							}
							else {
								Mix_FadeInChannelTimed(11, timerMusic, -1, 1500, -1);
								if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
							}
						}
						pauseSurface = false;
						t_continue.hover = false;
						if (firstFlag) startTime += (clock() - pauseTime);
					}
					else if (*x > t_back.x && *x<(t_back.x + t_back.w) && *y>t_back.y && *y < (t_back.y + t_back.h)) {
						Mix_PlayChannel(10, clickMusic, 0);
						Mix_HaltChannel(-1);
						menubgm = true;
						//Mix_PlayChannel(2, menuMusic, -1);
						t_back.hover = false;
						pauseSurface = false;
						firstFlag = false;
						firstDisplay = false;
						menuFlag = true;
						c_TIME.click = true;
						c_simple.click = true;
						//c_on.click = true;
						c_survival.click = false;
						c_INVERSION.click = false;
						c_classic.click = false;
						//c_off.click = false;
					}
				}
				*mouseState = IN_LB_SC;
			}
			else if (e->button.button == SDL_BUTTON_RIGHT && e->button.clicks == 1 && timediv < 800 && timediv > 30)
			{
				if (firstFlag && !(fail || succeed)) {
					userR -= M_PI;
				}

				*mouseState = IN_RB_SC;
			}

			break;

		case SDL_MOUSEWHEEL:
			if (e->wheel.y > 0) // scroll up
			{
				// Put code for handling "scroll up" here!
				*y = e->wheel.y;
				fFov--;
				*mouseState = IN_WU;
			}
			else if (e->wheel.y < 0) // scroll down
			{
				// Put code for handling "scroll down" here!
				*y = e->wheel.y;
				fFov++;
				*mouseState = IN_WD;
			}
			break;

		case SDL_MOUSEMOTION:

			*mouseState = HOVER;
			if (menuFlag) {

				if (*x > t_HOME.x && *x<(t_HOME.x + t_HOME.w) && *y>t_HOME.y && *y < (t_HOME.y + t_HOME.h)) {
					t_HOME.hover = true;
				}
				else t_HOME.hover = false;

				if (*x > c_TIME.x && *x<(c_TIME.x + c_TIME.w) && *y>c_TIME.y && *y < (c_TIME.y + c_TIME.h)) {
					c_TIME.hover = true;
				}
				else c_TIME.hover = false;

				if (*x > c_survival.x && *x<(c_survival.x + c_survival.w) && *y>c_survival.y && *y < (c_survival.y + c_survival.h)) {
					c_survival.hover = true;
				}
				else c_survival.hover = false;

				if (*x > c_simple.x && *x<(c_simple.x + c_simple.w) && *y>c_simple.y && *y < (c_simple.y + c_simple.h)) {
					c_simple.hover = true;
				}
				else c_simple.hover = false;

				if (*x > c_classic.x && *x<(c_classic.x + c_classic.w) && *y>c_classic.y && *y < (c_classic.y + c_classic.h)) {
					c_classic.hover = true;
				}
				else c_classic.hover = false;

				if (*x > c_on.x && *x<(c_on.x + c_on.w) && *y>c_on.y && *y < (c_on.y + c_on.h)) {
					c_on.hover = true;
				}
				else c_on.hover = false;

				if (*x > c_off.x && *x<(c_off.x + c_off.w) && *y>c_off.y && *y < (c_off.y + c_off.h)) {
					c_off.hover = true;
				}
				else c_off.hover = false;

				if (*x > t_PLAY.x && *x<(t_PLAY.x + t_PLAY.w) && *y>t_PLAY.y && *y < (t_PLAY.y + t_PLAY.h)) {
					t_PLAY.hover = true;
				}
				else t_PLAY.hover = false;

				if (*x > t_QUIT.x && *x<(t_QUIT.x + t_QUIT.w) && *y>t_QUIT.y && *y < (t_QUIT.y + t_QUIT.h)) {
					t_QUIT.hover = true;
				}
				else t_QUIT.hover = false;

				if (*x > c_INVERSION.x && *x<(c_INVERSION.x + c_INVERSION.w) && *y>c_INVERSION.y && *y < (c_INVERSION.y + c_INVERSION.h)) {
					c_INVERSION.hover = true;
				}
				else c_INVERSION.hover = false;
			}


			if ((firstFlag || jesusFlag) && !(fail || succeed) && sqrt(pow(*x - 100, 2) + powf((*y - 620), 2)) <= 50) {
				H_pause = true;
			}
			else H_pause = false;

			if ((fail || succeed) && *x > (t_AND.x + 12) && *x<(t_AND.x + (game_over_SizeW - 1) * 12) && *y>(t_AND.y + 12) && *y < (t_AND.y + 12 * (game_over_SizeH / 2 - 1))) {
				c_TryAgain.hover = true;
			}
			else c_TryAgain.hover = false;

			if ((fail || succeed) && *x > (t_AND.x + 12) && *x<(t_AND.x + (game_over_SizeW - 1) * 12) && *y>(t_AND.y + 12 * (game_over_SizeH / 2 + 1)) && *y < (t_AND.y + 12 * (game_over_SizeH - 1))) {
				c_BackToMenu.hover = true;
			}
			else c_BackToMenu.hover = false;

			if (pauseSurface) {
				if (*x > t_continue.x && *x<(t_continue.x + t_continue.w) && *y>t_continue.y && *y < (t_continue.y + t_continue.h)) {
					t_continue.hover = true;
				}
				else t_continue.hover = false;
				if (*x > t_back.x && *x<(t_back.x + t_back.w) && *y>t_back.y && *y < (t_back.y + t_back.h)) {
					t_back.hover = true;
				}
				else t_back.hover = false;
			}
			break;
		}

	}
}

void handleEvent(SDL_Event& e)
{
	if (e.type == SDL_KEYDOWN && e.key.repeat == 0)
	{
		static clock_t lastP = clock();
		static clock_t nowP = clock();
		int timediv;
		lastP = nowP;
		nowP = clock();
		timediv = nowP - lastP;
		//printf("timediv=%d\n", timediv);

		//first
		if (firstFlag && !pauseSurface) {
			/*if (userX < 0.5f || userX >= mapSize - 0.5f || userY < 0.5f || userY >= mapSize - 0.5f) {
				speed = 3.f;
			}
			else speed = 12.f;
			if (userX < 0.5f || userX >= mapSize - 0.5f || userY < 0.5f || userY >= mapSize - 0.5f) {
				fmod(userR, M_PI * 2.f);
				if (userR < 0) userR += 2 * M_PI;
				if (userR <= M_PI / 4.f && userR >= M_PI * 7.f / 4.f) userR = 0.f;
				else if (userR > M_PI / 4.f && userR < M_PI * 3.f / 4.f) userR = M_PI / 2.f;
				else if (userR >= M_PI * 3.f / 4.f && userR <= M_PI * 5.f / 4.f) userR = M_PI;
				else userR = M_PI * 3.f / 2.f;
			}*/

			switch (e.key.keysym.sym)
			{
			case SDLK_w:

				if (timediv > 50 && timediv < 500) {
					//printf("fast\n");
					shiftX = 2 * sinf(userR) * speed * elapsedTime;
					shiftY = 2 * cosf(userR) * speed * elapsedTime;
					//if (!quickFlag) quickFlag = true;
					//else quickFlag = false;
				}
				else {
					shiftX = sinf(userR) * speed * elapsedTime;
					shiftY = cosf(userR) * speed * elapsedTime;
				}
				break;

			case SDLK_s:
				if (timediv > 50 && timediv < 500) {
					//printf("fast\n");
					shiftX = -2 * sinf(userR) * speed * elapsedTime;
					shiftY = -2 * cosf(userR) * speed * elapsedTime;
					//if (!quickFlag) quickFlag = true;
					//else quickFlag = false;
				}
				else {
					shiftX = -sinf(userR) * speed * elapsedTime;
					shiftY = -cosf(userR) * speed * elapsedTime;
				}
				//shiftX = -sinf(userR) * speed * elapsedTime;
				//shiftY = -cosf(userR) * speed * elapsedTime;
				break;

			case SDLK_a:
				if (timediv > 50 && timediv < 500) {
					//printf("fast\n");
					shiftX = -2 * sinf(userR + M_PI / 2) * speed * elapsedTime;
					shiftY = -2 * cosf(userR + M_PI / 2) * speed * elapsedTime;
					//if (!quickFlag) quickFlag = true;
					//else quickFlag = false;
				}
				else {
					shiftX = -sinf(userR + M_PI / 2) * speed * elapsedTime;
					shiftY = -cosf(userR + M_PI / 2) * speed * elapsedTime;
				}
				//shiftX = -sinf(userR + M_PI / 2) * speed * elapsedTime;
				//shiftY = -cosf(userR + M_PI / 2) * speed * elapsedTime;
				break;

			case SDLK_d:
				if (timediv > 50 && timediv < 500) {
					//printf("fast\n");
					shiftX = 2 * sinf(userR + M_PI / 2) * speed * elapsedTime;
					shiftY = 2 * cosf(userR + M_PI / 2) * speed * elapsedTime;
					//if (!quickFlag) quickFlag = true;
					//else quickFlag = false;
				}
				else {
					shiftX = sinf(userR + M_PI / 2) * speed * elapsedTime;
					shiftY = cosf(userR + M_PI / 2) * speed * elapsedTime;
				}
				//shiftX = sinf(userR + M_PI / 2) * speed * elapsedTime;
				//shiftY = cosf(userR + M_PI / 2) * speed * elapsedTime;
				break;

			}
			/*if (quickFlag && !(userX < 1.f || userX >= mapSize - 1.f || userY < 1.f || userY >= mapSize - 1.f)) {
				shiftX *= 2.f;
				shiftY *= 2.f;
			}*/

		}

		switch (e.key.keysym.sym)
		{
		case SDLK_SPACE:
			Mix_PlayChannel(10, clickMusic, 0);
			if (entFlag && !menuFlag) {
				Mix_HaltChannel(1);
				//Mix_PlayChannel(2, menuMusic, -1);
				entFlag = false;
				menuFlag = true;
				menubgm = true;
				c_TIME.click = true;
				c_simple.click = true;
				c_on.click = true;
				c_survival.click = false;
				c_classic.click = false;
				c_off.click = false;
				c_INVERSION.click = false;
			}
			else if (transFlag) {
				transFlag = false;
				Mix_HaltChannel(9);
				if (c_survival.click) survivalbgm = true;
				else timerbgm = true;
				trans_start = 0;
				transFlag = false;
				firstFlag = true;
				firstDisplay = true;
				//jesusDisplay = true;
				duration = 3000; //change to 3 sec
			}
			break;
		case SDLK_TAB:
			//first
			if (firstFlag && !pauseSurface) {
				Mix_PlayChannel(10, clickMusic, 0);
				//Mix_HaltChannel(-1);
				if (c_survival.click) Mix_HaltChannel(5);
				else Mix_HaltChannel(11);
				Mix_PlayChannel(3, jesusMusic, -1);
				pauseTime = clock();
				firstFlag = false;
				jesusFlag = true;
				//jesusDisplay = true; //no need to initialize ?
			}

			//jesus
			else if (jesusFlag && !pauseSurface) {
				Mix_HaltChannel(3);
				if (c_survival.click) Mix_PlayChannel(5, survivalMusic, -1);
				else Mix_PlayChannel(11, timerMusic, -1);
				startTime += (clock() - pauseTime);
				firstFlag = true;
				jesusFlag = false;
				//firstDisplay = true;
			}
			break;


		case SDLK_p:
			//first 
			if (firstFlag && !(fail || succeed)) {
				if (!pauseSurface) {
					//if (c_survival.click) Mix_HaltChannel(5);
					//else Mix_HaltChannel(11);
					Mix_FadeOutChannel(-1, 1500);
					Mix_PlayChannel(6, pauseMusic, -1);
					pauseSurface = true;
					pauseTime = clock();
				}
				else {
					Mix_HaltChannel(-1);
					if (c_survival.click) {
						Mix_FadeInChannelTimed(5, survivalMusic, -1, 1500, -1);
						if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
					}
					else {
						Mix_FadeInChannelTimed(11, timerMusic, -1, 1500, -1);
						if (approach)Mix_FadeInChannelTimed(8, alertMusic, -1, 1500, -1);
					}
					pauseSurface = false;
					startTime += (clock() - pauseTime);
				}
			}
			//jesus
			else if (jesusFlag && !(fail || succeed)) {
				if (!pauseSurface) {
					Mix_HaltChannel(3);
					Mix_PlayChannel(6, pauseMusic, -1);
					pauseSurface = true;
				}
				else {
					pauseSurface = false;
					Mix_HaltChannel(6);
					Mix_PlayChannel(3, jesusMusic, -1);
				}
			}
			break;
		case SDLK_o:
			//jesus
			if (!firstFlag && jesusFlag && !pauseSurface) {
				recentFace = true;
			}
			break;
		}

	}
	/*
	else if (e.type == SDL_KEYDOWN && e.key.repeat == 1) {
		printf("in\n");
		static clock_t lastP = clock();
		static clock_t nowP = clock();
		int timediv;


		if (firstFlag && !pauseSurface) {
			switch (e.key.keysym.sym)
			{
			case SDLK_w:
				lastP = nowP;
				nowP = clock();
				timediv = nowP - lastP;
				if (timediv < 800) {
					if (!quickFlag) quickFlag = true;
					else quickFlag = false;

					shiftX *= 2;
					shiftY *= 2;

				}
				printf("time div = %d\n", timediv);
				break;

				case SDLK_s:
					shiftX *= 2;
					shiftY *= 2;

					break;

				case SDLK_a:
					shiftX *= 2;
					shiftY *= 2;
					break;
				case SDLK_d:
					shiftX *= 2;
					shiftY *= 2;
					break;

			}
		}
	}
	*/
	else if (e.type == SDL_KEYUP && e.key.repeat == 0)
	{
		switch (e.key.keysym.sym)
		{
		case SDLK_w:
			shiftX = 0;
			shiftY = 0;

			break;
		case SDLK_s:
			shiftX = 0;
			shiftY = 0;

			break;

		case SDLK_a:
			shiftX = 0;
			shiftY = 0;
			break;
		case SDLK_d:
			shiftX = 0;
			shiftY = 0;
			break;
		case SDLK_ESCAPE:
			if (menuFlag)
				exit(1);
		}

	}
}


