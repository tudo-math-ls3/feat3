// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>

namespace FEAT
{

/*
Copyright 2010 James Gregson. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY James Gregson ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL James Gregson OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/*
        A priority queue in which the priorities of entries may be
        updated in lg(n) time.

        Simply two maps glued together, one mapping keys (priorities)
        to values and the other mapping values to keys.  Both keys
        and values must meet the requirements for storage in a
        std::map, i.e. be less than comparable.  Typically not a
        problem for priorities, but can be a pain if Values are stored
        by value in the structures.  Not a problem if pointers to
        values are stored.

        insert, erase and update operations are all logarithmic, but
        approximately 2x slower than std::map since two std::maps
        are updated in each case.
*/
      template< typename _Key, typename _Value, typename _KeyComp=std::greater<_Key>, typename _ValueComp=std::less<_Value> >
      class mutable_priority_queue {
        public:
          // iterator type for the key2val map, allows in-order traversal of the queue
          typedef typename std::multimap<_Key,_Value,_KeyComp>::const_iterator iterator;
          typedef typename std::multimap<_Key,_Value,_KeyComp>::const_iterator const_iterator;

          // default constructor, uses std::greater and std::less for key and value comparisons
          mutable_priority_queue() : key2Val(_KeyComp()), val2Key(_ValueComp()) { }

          // constructor with key comparison predicate
          mutable_priority_queue( const _KeyComp &keyComp ) : key2Val(keyComp), val2Key(_ValueComp()){ }

          // constructor with value comparison predicate
          mutable_priority_queue( const _ValueComp &valComp ) : key2Val(_KeyComp()), val2Key(valComp) { }

          // constructor with both comparison predicates
          mutable_priority_queue( const _KeyComp &keyComp, const _ValueComp &valComp ) : key2Val(keyComp), val2Key(valComp) { }

          // copy constructor
          mutable_priority_queue( const mutable_priority_queue<_Key,_Value,_KeyComp,_ValueComp> &x ) : key2Val(x.key2Val), val2Key(x.val2Key) {}

          // destructor, clears both maps
          ~mutable_priority_queue(){
            key2Val.clear();
            val2Key.clear();
          }

          // assignment operator
          inline void operator=( mutable_priority_queue<_Key,_Value,_KeyComp,_ValueComp> x ){
            key2Val = x.key2Val;
            val2Key = x.val2Key;
          }

          // returns an iterator for the front of the queue
          inline iterator begin(){
            return key2Val.begin();
          }

          // returns an iterator for the end of the queue
          inline iterator end(){
            return key2Val.end();
          }

          // empties the queue
          inline void clear(){
            key2Val.clear();
            val2Key.clear();
          }

          // removes from both maps the value val
          inline void erase( _Value val ){
            typename std::map<_Value,_Key,_ValueComp>::iterator iter=val2Key.find( val );
            if( iter != val2Key.end() ){
              key2Val.erase( _find(iter->second,iter->first) );
              val2Key.erase( iter );
            }
          }

          // returns true if the queue is empty
          inline bool empty(){
            return key2Val.empty();
          }

          // returns the value at the front of the queue
          inline _Value front_value(){
            return (*key2Val.begin()).second;
          }

          // returns the key at the front of the queue
          inline _Value front_key(){
            return (*key2Val.begin()).first;
          }

          // removes the front entry in the queue
          inline void pop(){
            auto iter = key2Val.begin();
            _erase( iter );
          }

          // push an entry onto the queue
          inline void push( _Value val, _Key key ){
            insert( val, key );
          }

          // adds a value key pair to the queue.
          // an attempt first has to be made to erase
          // the item, since it could exist with another
          // key and the queue does not allow multiple
          // keys per value
          inline void insert( _Value val, _Key key ){
            erase( val );
            key2Val.insert( std::pair<_Key,_Value>( key, val ) );
            val2Key.insert( std::pair<_Value,_Key>( val, key ) );
          }

          // update the key associated with a value,
          // simply a convenience function that calls
          // insert()
          inline void update( _Value val, _Key key ){
            insert( val, key );
          }

          // returns the number of elements in the queue
          inline size_t size(){
            return key2Val.size();
          }

          inline size_t count(_Value val)
          {
            return val2Key.count(val);
          }

        private:
          std::multimap< _Key, _Value, _KeyComp > key2Val;
          std::map< _Value, _Key, _ValueComp >    val2Key;

          inline void _erase( typename std::multimap<_Key,_Value,_KeyComp>::iterator it ){
            val2Key.erase( it->second );
            key2Val.erase( it );
          }

          inline typename std::multimap<_Key,_Value,_KeyComp>::iterator _find( _Key k, _Value v ){
            typename std::multimap<_Key,_Value,_KeyComp>::iterator it = key2Val.lower_bound(k);
            typename std::multimap<_Key,_Value,_KeyComp>::iterator hi = key2Val.upper_bound(k);
            while( it != hi ){
              if( it->second == v )
                return it;
              it++;
            }
            return key2Val.end();
          }
      };
}
