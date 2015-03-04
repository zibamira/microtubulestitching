#include <hxalignmicrotubules/mtalign/Context.h>

#include <stdio.h>
#include <QString>

namespace ma = mtalign;

static void printStdout(QString msg) {
    puts(qPrintable(msg));
}

ma::Context::Context() : print(&printStdout) {}

ma::Context& ma::defaultContext() {
    static Context ctx;
    return ctx;
}
