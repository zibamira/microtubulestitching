#pragma once

#include <boost/function.hpp>

class QString;

namespace mtalign {

/// `Context` can be used to configure the run-time environment.
///
///  Several functions accept an optional `Context` object.  If provided,
///  the `Context` will control certain aspects of the function execution.
///
///  Set the function `print` to redirect printing.  `print` will be called
///  with a `QString` for each line (without terminating newline).
///
/// `Context` will probably be extended to also support `McProgressInterface`.
struct Context {
    Context();

    /// `print_t` is the type of a low-level print function.
    typedef boost::function<void (QString)> print_t;

    /// `print` is used to configure low-level printing.
    print_t print;
};

Context& defaultContext();

}  // namespace mtalign
