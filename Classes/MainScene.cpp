#include "MainScene.h"
#include "SimpleAudioEngine.h"
#include "TSPGA.h"

Scene* MainScene::createScene()
{
    return MainScene::create();
}

// on "init" you need to initialize your instance
bool MainScene::init()
{
    //////////////////////////////
    // 1. super init first
    if ( !Scene::init() )
    {
        return false;
    }

    auto visibleSize = Director::getInstance()->getVisibleSize();
    Vec2 origin = Director::getInstance()->getVisibleOrigin();

	auto keyListener = EventListenerKeyboard::create();
	keyListener->onKeyPressed = CC_CALLBACK_2(MainScene::OnKeyPressed, this);
	_eventDispatcher->addEventListenerWithSceneGraphPriority(keyListener, this);

	auto mouseListener = EventListenerMouse::create();
	mouseListener->onMouseDown = CC_CALLBACK_1(MainScene::OnMouseDown, this);
	_eventDispatcher->addEventListenerWithSceneGraphPriority(mouseListener, this);

	_drawNode = DrawNode::create();
	this->addChild(_drawNode);

	//show generation
	mp_labelGen = LabelTTF::create();
	mp_labelGen->setFontSize(16);
	mp_labelGen->setAnchorPoint(ccp(0, 0));
	mp_labelGen->setPosition(0, visibleSize.height - 20);
	mp_labelGen->setString(String::createWithFormat("Generation:%d",0)->getCString());
	this->addChild(mp_labelGen);

	//show information
	mp_labelInfo = LabelTTF::create();
	mp_labelInfo->setFontSize(16);
	mp_labelInfo->setAnchorPoint(ccp(0, 0.5));
	this->addChild(mp_labelInfo);
	mp_labelInfo->setPosition(0, mp_labelGen->getPositionY() - 30);
	tspga = new TSPGA(_cities, this);

	//choose crossover type button
	mp_menuItemLabel = LabelTTF::create();
	mp_menuItemLabel->setFontSize(16);
	mp_menuItemLabel->setString("crossover:PMX");
	auto menuItem = MenuItemLabel::create(mp_menuItemLabel, 
		CC_CALLBACK_1(MainScene::chooseCrossoverType, this));
	auto menu = Menu::create(menuItem, NULL);
	menu->setPosition(menuItem->getContentSize().width/2, 15);
	this->addChild(menu);

	scheduleUpdate();

	showInfo("-click Right Mouse to add city\n-click Left Mouse to move city");

    return true;
}


void MainScene::update(float dt) {
	if (tspga->mb_started && tspga->mb_running) {
		tspga->epoch();
		this->updateGeneration(tspga->mi_generation);
		this->clearPath();
		this->drawMap(tspga->getBestRoute());
	}
}

//add a city to the position where the right mouse clicked
void MainScene::OnMouseDown(Event *event) {
	EventMouse *me = (EventMouse*)event;
	auto btn = me->getMouseButton();
	if (btn == EventMouse::MouseButton::BUTTON_RIGHT) {
		//add a city
		auto city = City::create();
		//convert the OpenGL coordinate to cocos2d coordinate
		auto pos = me->getLocation();
		city->setPosition(pos.x, Director::getInstance()->getVisibleSize().height - pos.y);
		this->addChild(city);
		_cities.push_back(city);
	}
}


/*
* Enter:
	if you have set more than 15 cities and the genetic algorithm is not running,
	transfer the city array to the genetic algorithm to update it's TSPMap object
	and then start a new run.
  
  Esc:
	clear all the cities and the path between the cities, you have to reset the cities
	for a new run
  
  c/C:
	clear the path between the cities, you can start a new run without having to set 
	cities again. The new run will use the existing cities.
  
  Space:
	suspend/start the processing of genetic algorithm.
*/
void MainScene::OnKeyPressed(EventKeyboard::KeyCode code, Event* event) {
	if (code == EventKeyboard::KeyCode::KEY_ENTER) {
		//start a new run
		if (_cities.size() >= 15) {
			if (!tspga->mb_started) {
				tspga->setCities(_cities);
				tspga->createStartingPopulation();
				tspga->mb_started = true;
				tspga->mb_running = true;
				showInfo("-press SPACE to pause/continue\n-press ESC to cancel\n-press C to clear path");
			}
		}
		else {
			showInfo("number of cities should >= 15\n-press ESC to clear all cities\n-press ENTER to run");
		}
	}
	else if (code == EventKeyboard::KeyCode::KEY_ESCAPE) {
		//if (!tspga->mb_running) return;
		//clear all city
		for (int i = 0, n = _cities.size(); i < n; ++i) {
			_cities.at(i)->removeFromParentAndCleanup(true);
		}
		_cities.clear();
		clearPath();
		tspga->mb_started = false;
		showInfo("-click Right Mouse to add city\n-click Left Mouse to move city");
	}
	else if (code == EventKeyboard::KeyCode::KEY_SPACE) {
		//  pause / go
		if (!tspga->mb_started) return;
		tspga->mb_running = !tspga->mb_running;
		if (tspga->mb_running) tspga->mi_cntBestTime = 0;
		if (tspga->mb_running) 
			showInfo("-press SPACE to pause/continue\n-press ESC to clear all\n-press C to clear pathes");
	//	else 
	//		showInfo("-press SPACE to continue\n-press ESC to clear all\n-press C to clear pathes");
	}
	else if (code == EventKeyboard::KeyCode::KEY_C
		|| code == EventKeyboard::KeyCode::KEY_CAPITAL_C) {
		//clear path
		if (!tspga->mb_started) return;
		clearPath();
		tspga->mb_running = false;
		tspga->mb_started = false;
		showInfo("-press ENTER to start a new run\n-press ESC to clear all");
	}
}

//draw the city circle
void MainScene::drawMap(vector<Vec2> &cities) {
	int n = cities.size();
	for (int i = 0; i < n-1; ++i) {
		_drawNode->drawSegment( cities[i], cities[i+1], 1, Color4F::RED);
	}
	if(n>1) _drawNode->drawSegment(cities[n-1],cities[0],1,Color4F::RED);
}

void MainScene::updateGeneration(int n) {
	mp_labelGen->setString(String::createWithFormat("Generation:%d", n)->getCString());
}

void MainScene::showInfo(string info) {
	mp_labelInfo->setString(info);
}


void MainScene::chooseCrossoverType(Ref* ref) {
	int nextType = (tspga->getCrossoverType() + 1) % 4;
	switch (nextType) {
	case PMX: 
		mp_menuItemLabel->setString("crossover:PMX");
		break;
	case OBX: mp_menuItemLabel->setString("crossover:OBX");
		break;
	case PBX: mp_menuItemLabel->setString("crossover:PBX");
		break;
	case CroRandom: mp_menuItemLabel->setString("corssover:Random");
		break;
	}
	if (nextType == CroRandom) {
		tspga->mb_crossoverRandom = true;
	}
	else {		
		tspga->mb_crossoverRandom = false;
	}
	tspga->setCrossoverType(nextType);
}