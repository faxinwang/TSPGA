//#pragma once
#ifndef _CITY_
#define _CITY_
#include <cocos2d.h>
using namespace cocos2d;

class City :public Sprite {
public:
	CREATE_FUNC(City);
	bool init() {
		if (!Sprite::init()) return false;
		this->initWithFile("imgs/city.png");

		auto listener = EventListenerTouchOneByOne::create();
		listener->onTouchBegan = CC_CALLBACK_2(City::onTouchBegan, this);
		listener->onTouchMoved = CC_CALLBACK_2(City::onTouchMoved, this);
		_eventDispatcher->addEventListenerWithSceneGraphPriority(listener, this);

		return true;
	}

private:
	bool onTouchBegan(Touch* touch, Event*event) {
		return this->getBoundingBox().containsPoint(touch->getLocation());
	}
	void onTouchMoved(Touch* touch, Event*event) {
		auto visibleSize = Director::getInstance()->getVisibleSize();
		auto touchPos = touch->getLocation();
		if (touchPos.x < 0 || touchPos.x > visibleSize.width ||
			touchPos.y < 0 || touchPos.y > visibleSize.height) {
			return;
		}
		this->setPosition(touchPos);
	}
};


#endif