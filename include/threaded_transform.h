/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2015 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
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

#include <vector>
#include <cassert>

#include <pthread.h>

template <typename InputIterator>
struct thread_args {
  InputIterator first;
  InputIterator last;
};

template <typename Function, typename InputIterator>
void* trans(void* args_)
{
  typedef thread_args<InputIterator> argtype;

  argtype* args = reinterpret_cast<argtype*>(args_);

  InputIterator first = args->first, last = args->last;

  while (first != last) {
    *first = Function::op(*first);
    ++first;
  }

  pthread_exit(NULL);
}

template <typename Function, typename InputIterator>
void threaded_inplace(unsigned num_threads, InputIterator first,
                      InputIterator last) {
  using namespace std;

  typedef thread_args<InputIterator> argtype;

  int ret;
  unsigned i;

  size_t num_values = last - first;
  size_t num_values_per_threads = num_values / num_threads;

  vector<pthread_t> threads;
  vector<argtype> args;

  threads.resize(num_threads);
  args.resize(num_threads);

  for (i = 0; i < num_threads - 1; ++i) {
    args[i].first = first;
    args[i].last = first + num_values_per_threads;

    ret = pthread_create(&threads[i], NULL,
      trans<Function, InputIterator>, &args[i]);

    if (ret)
      break;

    first += num_values_per_threads;
  }

  if (ret) {
    assert(false);
  }

  // Last thread processes all remaining elements.
  args[i].first = first;
  args[i].last = last;
  ret = pthread_create(&threads[i], NULL,
    trans<Function, InputIterator>, &args[i]);

  for (vector<pthread_t>::iterator j = threads.begin(); j != threads.end(); ++j) {
    void *thread_ret;

    ret = pthread_join(*j, &thread_ret);

    if (ret)
      break;
  }

  if (ret) {
    assert(false);
  }
}
