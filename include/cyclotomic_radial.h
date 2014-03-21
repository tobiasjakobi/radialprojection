/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2014 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
 *
 *  radialprojection is free software: you can redistribute it and/or modify it under the terms
 *  of the GNU General Public License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  radialprojection is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with radialprojection.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CYCLOTOMIC_RADIAL_H_
#define _CYCLOTOMIC_RADIAL_H_

namespace SingleMachine {

  /* Default main routine for single machine execution */
  int main(int argc, char* argv[]);

};

namespace MultiMachine {

  /* multimachine routine that is executed on master machine */
  int master(int argc, char* argv[]);

  /* multimachine routine that is executed on the slave machines */
  int slave(int argc, char* argv[]);

};

#endif

