#ifndef __MainScene__
#define __MainScene__

#include "cocos2d.h"
#include <vector>

using namespace cocos2d;
using namespace std;

class City;
class TSPGA;

class MainScene : public cocos2d::Scene
{
public:
    static cocos2d::Scene* createScene();

    virtual bool init();
    
    // implement the "static create()" method manually
    CREATE_FUNC(MainScene);

	void OnMouseDown(Event* event);
	void OnKeyPressed(EventKeyboard::KeyCode code, Event* event);
	void drawMap(vector<Vec2> &cities);
	void clearPath() { _drawNode->clear(); }
	void updateGeneration(int n);
	void update(float dt);
	void showInfo(string info);

	void chooseCrossoverType(Ref* ref);

private:
	vector<City*> _cities;
	DrawNode *_drawNode;

	LabelTTF* mp_labelGen;
	LabelTTF* mp_labelInfo;
	LabelTTF* mp_menuItemLabel;

	TSPGA *tspga;

};

#endif // __HELLOWORLD_SCENE_H__
