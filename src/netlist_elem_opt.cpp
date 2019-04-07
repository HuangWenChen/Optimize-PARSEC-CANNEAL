// netlist_elem.cpp
//
// Created by Daniel Schwartz-Narbonne on 14/04/07.
//
// Copyright 2007 Princeton University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.


#include <iostream>
#include <assert.h>
#include <math.h>

#include "annealer_types.h"
#include "location_t.h"
#include "netlist_elem.h"

using namespace std;

netlist_elem::netlist_elem()
:present_loc(NULL)//start with the present_loc as nothing at all.  Filled in later by the netlist
{
}

//*****************************************************************************************
// Calculates the routing cost using the manhatten distance
// I make sure to get the pointer in one operation, and then use it
// SYNC: Do i need to make this an atomic operation?  i.e. are there misaligned memoery issues that can cause this to fail
//       even if I have atomic writes?
//*****************************************************************************************
routing_cost_t netlist_elem::routing_cost_given_loc(location_t loc)
{
	routing_cost_t fanin_cost = 0;
	routing_cost_t fanout_cost = 0;
	
	for (int i = 0; i< fanin.size(); ++i){
		location_t* fanin_loc = fanin[i]->present_loc.Get();
		fanin_cost += fabs(loc.x - fanin_loc->x);
		fanin_cost += fabs(loc.y - fanin_loc->y);
	}

	for (int i = 0; i< fanout.size(); ++i){
		location_t* fanout_loc = fanout[i]->present_loc.Get();
		fanout_cost += fabs(loc.x - fanout_loc->x);
		fanout_cost += fabs(loc.y - fanout_loc->y);
	}

	routing_cost_t total_cost = fanin_cost + fanout_cost;
	return total_cost;
}

//*****************************************************************************************
//  Get the cost change of swapping from our present location to a new location
//*****************************************************************************************
routing_cost_t netlist_elem::swap_cost(location_t* old_loc, location_t* new_loc)
{
	routing_cost_t no_swap = 0;
	routing_cost_t yes_swap = 0;
	
    long unsigned int remaining = fanin.size() % 4;
	int i = 0;
	
	for (i = 0; i< fanin.size() - remaining; i = i +4){
        location_t* fanin_loc1 = fanin[i]->present_loc.Get();
        location_t* fanin_loc2 = fanin[i+1]->present_loc.Get();
        location_t* fanin_loc3 = fanin[i+2]->present_loc.Get();
        location_t* fanin_loc4 = fanin[i+3]->present_loc.Get();
		no_swap += fabs(old_loc->x - fanin_loc1->x);
		no_swap += fabs(old_loc->y - fanin_loc1->y);
		no_swap += fabs(old_loc->x - fanin_loc2->x);
		no_swap += fabs(old_loc->y - fanin_loc2->y);
		no_swap += fabs(old_loc->x - fanin_loc3->x);
		no_swap += fabs(old_loc->y - fanin_loc3->y);
		no_swap += fabs(old_loc->x - fanin_loc4->x);
		no_swap += fabs(old_loc->y - fanin_loc4->y);
		
		yes_swap += fabs(new_loc->x - fanin_loc1->x);
		yes_swap += fabs(new_loc->y - fanin_loc1->y);
		yes_swap += fabs(new_loc->x - fanin_loc2->x);
		yes_swap += fabs(new_loc->y - fanin_loc2->y);
		yes_swap += fabs(new_loc->x - fanin_loc3->x);
		yes_swap += fabs(new_loc->y - fanin_loc3->y);
		yes_swap += fabs(new_loc->x - fanin_loc4->x);
		yes_swap += fabs(new_loc->y - fanin_loc4->y);
	}
	
    if(remaining){
		for (; i< fanin.size(); i++){
			location_t* fanin_loc1 = fanin[i]->present_loc.Get();
			no_swap += fabs(old_loc->x - fanin_loc1->x);
			no_swap += fabs(old_loc->y - fanin_loc1->y);
			
			yes_swap += fabs(new_loc->x - fanin_loc1->x);
			yes_swap += fabs(new_loc->y - fanin_loc1->y);
		}
    }


    remaining = fanout.size() % 4;
	i = 0;
	
	for (i = 0; i< fanout.size() - remaining; i = i +4){
        location_t* fanout_loc1 = fanout[i]->present_loc.Get();
        location_t* fanout_loc2 = fanout[i+1]->present_loc.Get();
        location_t* fanout_loc3 = fanout[i+2]->present_loc.Get();
        location_t* fanout_loc4 = fanout[i+3]->present_loc.Get();
		no_swap += fabs(old_loc->x - fanout_loc1->x);
		no_swap += fabs(old_loc->y - fanout_loc1->y);
		no_swap += fabs(old_loc->x - fanout_loc2->x);
		no_swap += fabs(old_loc->y - fanout_loc2->y);
		no_swap += fabs(old_loc->x - fanout_loc3->x);
		no_swap += fabs(old_loc->y - fanout_loc3->y);
		no_swap += fabs(old_loc->x - fanout_loc4->x);
		no_swap += fabs(old_loc->y - fanout_loc4->y);
		
		yes_swap += fabs(new_loc->x - fanout_loc1->x);
		yes_swap += fabs(new_loc->y - fanout_loc1->y);
		yes_swap += fabs(new_loc->x - fanout_loc2->x);
		yes_swap += fabs(new_loc->y - fanout_loc2->y);
		yes_swap += fabs(new_loc->x - fanout_loc3->x);
		yes_swap += fabs(new_loc->y - fanout_loc3->y);
		yes_swap += fabs(new_loc->x - fanout_loc4->x);
		yes_swap += fabs(new_loc->y - fanout_loc4->y);
	}
	
    if(remaining){
		for (; i< fanout.size(); i++){
			location_t* fanout_loc1 = fanout[i]->present_loc.Get();
			no_swap += fabs(old_loc->x - fanout_loc1->x);
			no_swap += fabs(old_loc->y - fanout_loc1->y);
			
			yes_swap += fabs(new_loc->x - fanout_loc1->x);
			yes_swap += fabs(new_loc->y - fanout_loc1->y);
		}
    }
	
	return yes_swap - no_swap;
}

